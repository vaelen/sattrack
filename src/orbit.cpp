/*
 * Copyright (c) 2025 Andrew C. Young <andrew@vaelen.org>
 * SPDX-License-Identifier: MIT
 */

#include <sattrack/orbit.hpp>

#include <charconv>
#include <chrono>
#include <cmath>
#include <numbers>
#include <ranges>
#include <string>
#include <string_view>
#include <filesystem>
#include <fstream>

#include <date/date.h>

namespace sattrack {

// Helper function to trim leading spaces from a string_view
std::string_view trimLeft(const std::string_view &str) {
    auto pos = str.find_first_not_of(' ');
    return pos == std::string_view::npos ? "" : str.substr(pos);
}

// Helper function to trim trailing spaces from a string_view
std::string_view trimRight(const std::string_view &str) {
    auto pos = str.find_last_not_of(' ');
    return pos == std::string_view::npos ? "" : str.substr(0, pos + 1);
}

// Helper function to convert substring to numeric type
template <typename T>
inline T toNumber(const std::string_view &str) {
    T value;
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), value);
    if (ec != std::errc()) {
        throw std::invalid_argument("Couldn't convert value: " + std::string(str));
    }
    return value;
}

// Helper function to convert exponential substring to numeric type
inline double fromExponentialString(const std::string_view &str) {
    // Example input: "11606-4" -> 0.00011606

    bool isNegative = true;
    auto pos = str.find('-');
    if (pos == std::string_view::npos) {
        isNegative = false;
        pos = str.find('+');
        if (pos == std::string_view::npos) {
            throw std::invalid_argument("Invalid exponential format: " + std::string(str));
        }
    }

    std::string baseStr = "0." + std::string(str.substr(0, pos));
    std::string_view exponentView = str.substr(pos + 1);

    double base = toNumber<double>(baseStr);
    int exponent = toNumber<int>(exponentView);
    double value = base * std::pow(10.0, isNegative ? -exponent : exponent);
    return value;
}

// Helper function to parse epoch from TLE format (YYDDD.DDDDDDDD)
auto parseEpoch(const std::string_view &epochStr) {
    using namespace std::chrono;

    int y = toNumber<int>(trimLeft(epochStr.substr(0, 2)));
    double dayOfYear = toNumber<double>(trimLeft(epochStr.substr(2)));

    // Convert two-digit year to four-digit year
    if (y < 57) {
        y += 2000;
    } else {
        y += 1900;
    }

    int wholeDays = static_cast<int>(dayOfYear);
    double fracDays = dayOfYear - wholeDays;
    
    auto date = sys_days{year{y}/January/1} + days{wholeDays - 1};
    auto time = duration_cast<microseconds>(duration<double, std::ratio<86400>>{fracDays});
    
    return date + time;
}

// Convert a time_point to Julian Date
double toJulianDate(time_point tp) {
    using namespace std::chrono;
    
    auto daysSinceEpoch = duration_cast<duration<double, days::period>>(
        tp.time_since_epoch()
    ).count();
    
    // Unix epoch (1970-01-01) in Julian Date is 2440587.5
    return 2440587.5 + daysSinceEpoch;
}

// Greenwich Mean Sidereal Time in radians
double gmst(double julianDate) {
    // Julian centuries since J2000.0
    double T = (julianDate - J2000_JD) / DAYS_PER_JULIAN_CENTURY;

    // GMST in degrees at 0h UT, then add rotation for time of day
    double gmstInDegrees = GMST_AT_J2000
                    + EARTH_SIDEREAL_RATE * (julianDate - J2000_JD)
                    + GMST_T2_COEFF * T * T
                    - T * T * T / GMST_T3_DIVISOR;
    
    // Normalize to [0, 360)
    gmstInDegrees = std::fmod(gmstInDegrees, 360.0);
    if (gmstInDegrees < 0) gmstInDegrees += 360.0;
    
    return gmstInDegrees * DEGREES_TO_RADIANS;
}

// ============================================================================
// SGP4/SDP4 Implementation
// Based on the Vallado reference implementation from CelesTrak
// ============================================================================

// Helper: Compute Greenwich Sidereal Time at epoch for SGP4
static double gstime(double jdut1) {
    double tut1 = (jdut1 - 2451545.0) / 36525.0;
    double temp = -6.2e-6 * tut1 * tut1 * tut1
                  + 0.093104 * tut1 * tut1
                  + (876600.0 * 3600 + 8640184.812866) * tut1
                  + 67310.54841;
    temp = std::fmod(temp * DEGREES_TO_RADIANS / 240.0, SGP4_TWO_PI);
    if (temp < 0.0) temp += SGP4_TWO_PI;
    return temp;
}

// Initialize SGP4 state from orbital elements
void Orbit::initializeSGP4() const {
    // Convert TLE elements to SGP4 internal units
    ecco_ = eccentricity;
    inclo_ = inclination * DEGREES_TO_RADIANS;
    nodeo_ = rightAscensionOfAscendingNode * DEGREES_TO_RADIANS;
    argpo_ = argumentOfPerigee * DEGREES_TO_RADIANS;
    mo_ = meanAnomaly * DEGREES_TO_RADIANS;
    bstar_ = bstarDragTerm;

    // Convert mean motion from revs/day to rad/min
    no_kozai_ = meanMotion * SGP4_TWO_PI / 1440.0;

    // Compute epoch Julian Date
    double epochJD = toJulianDate(epoch);
    jdsatepoch_ = std::floor(epochJD);
    jdsatepochF_ = epochJD - jdsatepoch_;

    // Compute Greenwich sidereal time at epoch
    gsto_ = gstime(epochJD);

    // WGS-72 Earth constants
    constexpr double radiusearthkm = SGP4_RADIUS_EARTH_KM;
    constexpr double xke = SGP4_XKE;
    constexpr double j2 = SGP4_J2;
    constexpr double j3oj2 = SGP4_J3OJ2;
    constexpr double j4 = SGP4_J4;

    // Recover original mean motion (no_unkozai) and semimajor axis from input
    double a1 = std::pow(xke / no_kozai_, SGP4_X2O3);
    double cosio = std::cos(inclo_);
    double cosio2 = cosio * cosio;
    double eccsq = ecco_ * ecco_;
    double omeosq = 1.0 - eccsq;
    double rteosq = std::sqrt(omeosq);
    double d1 = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
    double del_ = d1 / (a1 * a1);
    double ao = a1 * (1.0 - del_ * (1.0/3.0 + del_ * (1.0 + 134.0/81.0 * del_)));
    double delo = d1 / (ao * ao);
    no_unkozai_ = no_kozai_ / (1.0 + delo);
    ao = a1 * (1.0 - del_ * (1.0/3.0 + del_ * (1.0 + 134.0/81.0 * del_)));

    // Compute semi-major axis
    a_ = std::pow(xke / no_unkozai_, SGP4_X2O3);

    // Compute perigee and apogee altitudes
    altp_ = a_ * (1.0 - ecco_) - 1.0;
    alta_ = a_ * (1.0 + ecco_) - 1.0;

    // Determine if deep space (period >= 225 min)
    double periodearthradii = SGP4_TWO_PI / no_unkozai_;
    if (periodearthradii >= 225.0) {
        method_ = 'd';  // deep space
    } else {
        method_ = 'n';  // near earth
    }

    // SGP4 initialization
    double ss = 78.0 / radiusearthkm + 1.0;
    double qzms2t = std::pow((120.0 - 78.0) / radiusearthkm, 4);

    double sinio = std::sin(inclo_);
    double x1mth2_ = 1.0 - cosio2;
    this->x1mth2_ = x1mth2_;
    cosio2_ = cosio2;

    // Check for eccentricity out of range
    if (ecco_ >= 1.0 || ecco_ < -0.001) {
        throw InvalidOrbitException("Eccentricity out of range: " + std::to_string(ecco_));
    }
    if (ecco_ < 1.0e-10) {
        ecco_ = 1.0e-10;  // Clamp to avoid division by zero
    }

    // Compute perigee
    double rp = a_ * (1.0 - ecco_);

    // Check if satellite has decayed
    if (rp < 1.0) {
        throw SatelliteDecayedException();
    }

    double sfour = ss;
    double qzms24 = qzms2t;
    double perige = (rp - 1.0) * radiusearthkm;

    // For perigees below 156 km, adjust s and qoms2t
    if (perige < 156.0) {
        sfour = perige - 78.0;
        if (perige < 98.0) {
            sfour = 20.0;
        }
        qzms24 = std::pow((120.0 - sfour) / radiusearthkm, 4);
        sfour = sfour / radiusearthkm + 1.0;
    }

    double pinvsq = 1.0 / (a_ * a_ * omeosq * omeosq);
    double tsi = 1.0 / (a_ - sfour);
    eta_ = a_ * ecco_ * tsi;
    double etasq = eta_ * eta_;
    double eeta = ecco_ * eta_;
    double psisq = std::abs(1.0 - etasq);
    double coef = qzms24 * std::pow(tsi, 4);
    double coef1 = coef / std::pow(psisq, 3.5);
    double cc2 = coef1 * no_unkozai_ * (a_ * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
                 + 0.375 * j2 * tsi / psisq * (3.0 * (3.0 * cosio2 - 1.0)
                 * (1.0 + 1.5 * etasq) - 0.5 * (3.0 - 7.0 * cosio2) * eeta * std::cos(2.0 * argpo_)));
    cc1_ = bstar_ * cc2;
    double cc3 = 0.0;
    if (ecco_ > 1.0e-4) {
        cc3 = -2.0 * coef * tsi * j3oj2 * no_unkozai_ * sinio / ecco_;
    }
    x7thm1_ = 7.0 * cosio2 - 1.0;
    cc4_ = 2.0 * no_unkozai_ * coef1 * a_ * omeosq
           * (eta_ * (2.0 + 0.5 * etasq) + ecco_ * (0.5 + 2.0 * etasq)
           - j2 * tsi / (a_ * psisq) * (-3.0 * (3.0 * (3.0 * cosio2 - 1.0) - 7.0 * x7thm1_)
           * (1.0 + 1.5 * etasq) / 6.0 + 0.25 * (3.0 - 7.0 * cosio2) * (2.0 * etasq - eeta * (1.0 + etasq)) * std::cos(2.0 * argpo_)));
    cc5_ = 2.0 * coef1 * a_ * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);

    double cosio4 = cosio2 * cosio2;
    double temp1 = 1.5 * j2 * pinvsq * no_unkozai_;
    double temp2 = 0.5 * temp1 * j2 * pinvsq;
    double temp3 = -0.46875 * j4 * pinvsq * pinvsq * no_unkozai_;
    mdot_ = no_unkozai_ + 0.5 * temp1 * rteosq * (3.0 * cosio2 - 1.0)
            + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
    argpdot_ = -0.5 * temp1 * (1.0 - 5.0 * cosio2)
               + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4)
               + temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
    double xhdot1 = -temp1 * cosio;
    nodedot_ = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
    omgcof_ = bstar_ * cc3 * std::cos(argpo_);
    xmcof_ = 0.0;
    if (ecco_ > 1.0e-4) {
        xmcof_ = -SGP4_X2O3 * coef * bstar_ / eeta;
    }
    nodecf_ = 3.5 * omeosq * xhdot1 * cc1_;
    t2cof_ = 1.5 * cc1_;

    // Set xlcof and aycof
    if (std::abs(cosio + 1.0) > 1.5e-12) {
        xlcof_ = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    } else {
        xlcof_ = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / 1.5e-12;
    }
    aycof_ = -0.5 * j3oj2 * sinio;

    // For SGP4, initialize additional terms
    delmo_ = std::pow(1.0 + eta_ * std::cos(mo_), 3);
    sinmao_ = std::sin(mo_);

    // Compute con41 for later use
    con41_ = 3.0 * cosio2 - 1.0;

    // Set isimp flag for very high drag satellites
    isimp_ = false;
    if ((omeosq >= 0.0) || (no_unkozai_ >= 0.0)) {
        if ((rp < (220.0 / radiusearthkm + 1.0))) {
            isimp_ = true;
        }
    }

    // Initialize d2, d3, d4, t3cof, t4cof, t5cof (for non-simple satellites)
    d2_ = 0.0;
    d3_ = 0.0;
    d4_ = 0.0;
    t3cof_ = 0.0;
    t4cof_ = 0.0;
    t5cof_ = 0.0;

    if (!isimp_) {
        double c1sq = cc1_ * cc1_;
        d2_ = 4.0 * a_ * tsi * c1sq;
        double temp = d2_ * tsi * cc1_ / 3.0;
        d3_ = (17.0 * a_ + sfour) * temp;
        d4_ = 0.5 * temp * a_ * tsi * (221.0 * a_ + 31.0 * sfour) * cc1_;
        t3cof_ = d2_ + 2.0 * c1sq;
        t4cof_ = 0.25 * (3.0 * d3_ + cc1_ * (12.0 * d2_ + 10.0 * c1sq));
        t5cof_ = 0.2 * (3.0 * d4_ + 12.0 * cc1_ * d3_ + 6.0 * d2_ * d2_ + 15.0 * c1sq * (2.0 * d2_ + c1sq));
    }

    // Deep space initialization if needed
    if (method_ == 'd') {
        irez_ = 0;
        if ((no_unkozai_ < 0.0052359877) && (no_unkozai_ > 0.0034906585)) {
            irez_ = 1;  // Synchronous resonance
        }
        if ((no_unkozai_ >= 8.26e-3) && (no_unkozai_ <= 9.24e-3) && (ecco_ >= 0.5)) {
            irez_ = 2;  // Half-day resonance
        }

        // Initialize deep space terms
        double tc = 0.0;
        double snodm, cnodm, sinim, cosim, sinomm, cosomm;
        double day, em, emsq, gam, rtemsq;
        double s1, s2, s3, s4, s5, s6, s7;
        double ss1, ss2, ss3, ss4, ss5, ss6, ss7;
        double sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33;
        double nm, z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33;

        initializeDeepSpace(tc, snodm, cnodm, sinim, cosim, sinomm, cosomm,
                           day, em, emsq, gam, rtemsq, s1, s2, s3, s4, s5, s6, s7,
                           ss1, ss2, ss3, ss4, ss5, ss6, ss7,
                           sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33,
                           nm, z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33);
    }

    sgp4Initialized_ = true;
}

// Initialize deep space (SDP4) coefficients
void Orbit::initializeDeepSpace(
    double tc, double& snodm, double& cnodm, double& sinim, double& cosim,
    double& sinomm, double& cosomm, double& day, double& em, double& emsq,
    double& gam, double& rtemsq, double& s1, double& s2, double& s3,
    double& s4, double& s5, double& s6, double& s7, double& ss1,
    double& ss2, double& ss3, double& ss4, double& ss5, double& ss6,
    double& ss7, double& sz1, double& sz2, double& sz3, double& sz11,
    double& sz12, double& sz13, double& sz21, double& sz22, double& sz23,
    double& sz31, double& sz32, double& sz33, double& nm, double& z1,
    double& z2, double& z3, double& z11, double& z12, double& z13,
    double& z21, double& z22, double& z23, double& z31, double& z32,
    double& z33) const {

    // Solar and lunar constants
    constexpr double zes = 0.01675;
    constexpr double zel = 0.05490;
    constexpr double c1ss = 2.9864797e-6;
    constexpr double c1l = 4.7968065e-7;
    constexpr double zsinis = 0.39785416;
    constexpr double zcosis = 0.91744867;
    constexpr double zcosgs = 0.1945905;
    constexpr double zsings = -0.98088458;

    nm = no_unkozai_;
    em = ecco_;
    snodm = std::sin(nodeo_);
    cnodm = std::cos(nodeo_);
    sinomm = std::sin(argpo_);
    cosomm = std::cos(argpo_);
    sinim = std::sin(inclo_);
    cosim = std::cos(inclo_);
    emsq = em * em;
    double betasq = 1.0 - emsq;
    rtemsq = std::sqrt(betasq);

    day = jdsatepoch_ + jdsatepochF_ - 2433281.5 + tc / 1440.0;
    double xnodce = std::fmod(4.5236020 - 9.2422029e-4 * day, SGP4_TWO_PI);
    double stem = std::sin(xnodce);
    double ctem = std::cos(xnodce);
    double zcosil = 0.91375164 - 0.03568096 * ctem;
    double zsinil = std::sqrt(1.0 - zcosil * zcosil);
    double zsinhl = 0.089683511 * stem / zsinil;
    double zcoshl = std::sqrt(1.0 - zsinhl * zsinhl);
    gam = 5.8351514 + 0.0019443680 * day;
    double zx = 0.39785416 * stem / zsinil;
    double zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    zx = std::atan2(zx, zy);
    zx = gam + zx - xnodce;
    double zcosgl = std::cos(zx);
    double zsingl = std::sin(zx);

    // Solar terms
    double zcosg = zcosgs;
    double zsing = zsings;
    double zcosi = zcosis;
    double zsini = zsinis;
    double zcosh = cnodm;
    double zsinh = snodm;
    double cc = c1ss;
    double xnoi = 1.0 / nm;

    for (int lsflg = 1; lsflg <= 2; lsflg++) {
        double a1 = zcosg * zcosh + zsing * zcosi * zsinh;
        double a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
        double a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
        double a8 = zsing * zsini;
        double a9 = zsing * zsinh + zcosg * zcosi * zcosh;
        double a10 = zcosg * zsini;
        double a2 = cosim * a7 + sinim * a8;
        double a4 = cosim * a9 + sinim * a10;
        double a5 = -sinim * a7 + cosim * a8;
        double a6 = -sinim * a9 + cosim * a10;

        double x1 = a1 * cosomm + a2 * sinomm;
        double x2 = a3 * cosomm + a4 * sinomm;
        double x3 = -a1 * sinomm + a2 * cosomm;
        double x4 = -a3 * sinomm + a4 * cosomm;
        double x5 = a5 * sinomm;
        double x6 = a6 * sinomm;
        double x7 = a5 * cosomm;
        double x8 = a6 * cosomm;

        z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
        z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
        z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
        z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * emsq;
        z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * emsq;
        z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * emsq;
        z11 = -6.0 * a1 * a5 + emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
        z12 = -6.0 * (a1 * a6 + a3 * a5) + emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
        z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
        z21 = 6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
        z22 = 6.0 * (a4 * a5 + a2 * a6) + emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
        z23 = 6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
        z1 = z1 + z1 + betasq * z31;
        z2 = z2 + z2 + betasq * z32;
        z3 = z3 + z3 + betasq * z33;
        s3 = cc * xnoi;
        s2 = -0.5 * s3 / rtemsq;
        s4 = s3 * rtemsq;
        s1 = -15.0 * em * s4;
        s5 = x1 * x3 + x2 * x4;
        s6 = x2 * x3 + x1 * x4;
        s7 = x2 * x4 - x1 * x3;

        if (lsflg == 1) {
            ss1 = s1;
            ss2 = s2;
            ss3 = s3;
            ss4 = s4;
            ss5 = s5;
            ss6 = s6;
            ss7 = s7;
            sz1 = z1;
            sz2 = z2;
            sz3 = z3;
            sz11 = z11;
            sz12 = z12;
            sz13 = z13;
            sz21 = z21;
            sz22 = z22;
            sz23 = z23;
            sz31 = z31;
            sz32 = z32;
            sz33 = z33;
            zcosg = zcosgl;
            zsing = zsingl;
            zcosi = zcosil;
            zsini = zsinil;
            zcosh = zcoshl * cnodm + zsinhl * snodm;
            zsinh = snodm * zcoshl - cnodm * zsinhl;
            cc = c1l;
        }
    }

    zmol_ = std::fmod(4.7199672 + 0.22997150 * day - gam, SGP4_TWO_PI);
    zmos_ = std::fmod(6.2565837 + 0.017201977 * day, SGP4_TWO_PI);

    // Solar terms
    se2_ = 2.0 * ss1 * ss6;
    se3_ = 2.0 * ss1 * ss7;
    si2_ = 2.0 * ss2 * sz12;
    si3_ = 2.0 * ss2 * (sz13 - sz11);
    sl2_ = -2.0 * ss3 * sz2;
    sl3_ = -2.0 * ss3 * (sz3 - sz1);
    sl4_ = -2.0 * ss3 * (-21.0 - 9.0 * emsq) * zes;
    sgh2_ = 2.0 * ss4 * sz32;
    sgh3_ = 2.0 * ss4 * (sz33 - sz31);
    sgh4_ = -18.0 * ss4 * zes;
    sh2_ = -2.0 * ss2 * sz22;
    sh3_ = -2.0 * ss2 * (sz23 - sz21);

    // Lunar terms
    ee2_ = 2.0 * s1 * s6;
    e3_ = 2.0 * s1 * s7;
    xi2_ = 2.0 * s2 * z12;
    xi3_ = 2.0 * s2 * (z13 - z11);
    xl2_ = -2.0 * s3 * z2;
    xl3_ = -2.0 * s3 * (z3 - z1);
    xl4_ = -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel;
    xgh2_ = 2.0 * s4 * z32;
    xgh3_ = 2.0 * s4 * (z33 - z31);
    xgh4_ = -18.0 * s4 * zel;
    xh2_ = -2.0 * s2 * z22;
    xh3_ = -2.0 * s2 * (z23 - z21);

    // Apply deep space long period periodics
    double f2 = 0.5 * sinomm * sinomm - 0.25;
    double f3 = -0.5 * sinomm * cosomm;
    double ses = se2_ * f2 + se3_ * f3;
    double sis = si2_ * f2 + si3_ * f3;
    double sls = sl2_ * f2 + sl3_ * f3 + sl4_ * sinomm;
    double sghs = sgh2_ * f2 + sgh3_ * f3 + sgh4_ * sinomm;
    double shs = sh2_ * f2 + sh3_ * f3;

    if (inclo_ < 5.2359877e-2 || inclo_ > M_PI - 5.2359877e-2) {
        shs = 0.0;
    }
    if (sinim != 0.0) {
        shs = shs / sinim;
    }

    double sel = ee2_ * f2 + e3_ * f3;
    double sil = xi2_ * f2 + xi3_ * f3;
    double sll = xl2_ * f2 + xl3_ * f3 + xl4_ * sinomm;
    double sghl = xgh2_ * f2 + xgh3_ * f3 + xgh4_ * sinomm;
    double shll = xh2_ * f2 + xh3_ * f3;

    if (inclo_ < 5.2359877e-2 || inclo_ > M_PI - 5.2359877e-2) {
        shll = 0.0;
    }

    peo_ = ses + sel;
    pinco_ = sis + sil;
    plo_ = sls + sll;
    pgho_ = sghs + sghl;
    pho_ = shs + shll;

    // Initialize resonance terms
    if (irez_ != 0) {
        double aonv = std::pow(nm / SGP4_XKE, SGP4_X2O3);

        // Compute the resonance terms
        // Geopotential resonance terms
        double cosisq = cosim * cosim;
        double emo = em;
        em = ecco_;
        double emsqo = emsq;
        emsq = eccsq_;
        emsq = ecco_ * ecco_;
        double eoc = em * emsq;
        double g201 = -0.306 - (em - 0.64) * 0.440;
        double g211, g310, g322, g410, g422, g520, g521, g532, g533;

        if (em <= 0.65) {
            g211 = 3.616 - 13.2470 * em + 16.2900 * emsq;
            g310 = -19.302 + 117.3900 * em - 228.4190 * emsq + 156.5910 * eoc;
            g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq + 146.5816 * eoc;
            g410 = -41.122 + 242.6940 * em - 471.0940 * emsq + 313.9530 * eoc;
            g422 = -146.407 + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
            g520 = -532.114 + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
        } else {
            g211 = -72.099 + 331.819 * em - 508.738 * emsq + 266.724 * eoc;
            g310 = -346.844 + 1582.851 * em - 2415.925 * emsq + 1246.113 * eoc;
            g322 = -342.585 + 1554.908 * em - 2366.899 * emsq + 1215.972 * eoc;
            g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq + 3651.957 * eoc;
            g422 = -3581.690 + 16178.11 * em - 24462.77 * emsq + 12422.52 * eoc;
            if (em > 0.715) {
                g520 = -5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
            } else {
                g520 = 1464.74 - 4664.75 * em + 3763.64 * emsq;
            }
        }

        if (em < 0.7) {
            g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.210 * eoc;
            g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
            g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.400 * eoc;
        } else {
            g533 = -37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
            g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
            g532 = -40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
        }

        double sini2 = sinim * sinim;
        double f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq);
        double f221 = 1.5 * sini2;
        double f321 = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
        double f322 = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
        double f441 = 35.0 * sini2 * f220;
        double f442 = 39.375 * sini2 * sini2;
        double f522 = 9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq) +
                      0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
        double f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim + 10.0 * cosisq) +
                      6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
        double f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
        double f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));

        double xno2 = nm * nm;
        double ainv2 = aonv * aonv;
        double temp1_ds = 3.0 * xno2 * ainv2;
        double temp_ds = temp1_ds * SGP4_J2 / 2.0;

        if (irez_ == 2) {
            // Half-day resonance terms
            double theta = std::fmod(gsto_ + tc * 7.29211514668855e-5, SGP4_TWO_PI);
            d2201_ = temp_ds * f220 * g201;
            d2211_ = temp_ds * f221 * g211;
            d3210_ = temp1_ds * 1.5 * SGP4_J2 * f321 * g310 * ainv2;
            d3222_ = temp1_ds * 1.5 * SGP4_J2 * f322 * g322 * ainv2;
            d4410_ = temp1_ds * temp_ds * 2.0 * f441 * g410 * ainv2;
            d4422_ = temp1_ds * temp_ds * 2.0 * f442 * g422 * ainv2;
            d5220_ = temp1_ds * temp_ds * f522 * g520 * ainv2 * ainv2;
            d5232_ = temp1_ds * temp_ds * f523 * g532 * ainv2 * ainv2;
            d5421_ = temp1_ds * temp_ds * f542 * g521 * ainv2 * ainv2;
            d5433_ = temp1_ds * temp_ds * f543 * g533 * ainv2 * ainv2;
            xlamo_ = std::fmod(mo_ + nodeo_ + nodeo_ - theta - theta, SGP4_TWO_PI);
            double xfact = mdot_ + dmdt_ + 2.0 * (nodedot_ + dnodt_ - 7.29211514668855e-5) - no_unkozai_;
            em = emo;
            emsq = emsqo;
        }

        if (irez_ == 1) {
            // Synchronous resonance terms
            double g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
            double g310_s = 1.0 + 2.0 * emsq;
            double g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
            double f220_s = 0.75 * (1.0 + cosim) * (1.0 + cosim);
            double f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
            double f330 = 1.0 + cosim;
            f330 = 1.875 * f330 * f330 * f330;

            del1_ = 3.0 * nm * nm * aonv * aonv;
            del2_ = 2.0 * del1_ * f220_s * g200 * SGP4_J2;
            del3_ = 3.0 * del1_ * f330 * g300 * SGP4_J3OJ2 * aonv;
            del1_ = del1_ * f311 * g310_s * SGP4_J2 * aonv;

            double theta_s = std::fmod(gsto_ + tc * 7.29211514668855e-5, SGP4_TWO_PI);
            xlamo_ = std::fmod(mo_ + nodeo_ + argpo_ - theta_s, SGP4_TWO_PI);
            double xfact = mdot_ + argpdot_ + nodedot_ - 7.29211514668855e-5 - no_unkozai_;
        }

        // Initialize integrator
        xli_ = xlamo_;
        xni_ = no_unkozai_;
        atime_ = 0.0;
    }
}

// Apply deep space secular effects
void Orbit::deepSpaceSecular(
    double t, double& em, double& argpm, double& inclm,
    double& nodem, double& mm, double& nm) const {

    constexpr double step = 720.0;
    constexpr double step2 = step * step / 2.0;

    // Apply lunar-solar periodics
    double zm = zmos_ + 0.017201977 * t;
    double zf = zm + 2.0 * 0.01675 * std::sin(zm);
    double sinzf = std::sin(zf);
    double f2 = 0.5 * sinzf * sinzf - 0.25;
    double f3 = -0.5 * sinzf * std::cos(zf);

    double ses = se2_ * f2 + se3_ * f3;
    double sis = si2_ * f2 + si3_ * f3;
    double sls = sl2_ * f2 + sl3_ * f3 + sl4_ * sinzf;
    double sghs = sgh2_ * f2 + sgh3_ * f3 + sgh4_ * sinzf;
    double shs = sh2_ * f2 + sh3_ * f3;

    zm = zmol_ + 0.22997150 * t;
    zf = zm + 2.0 * 0.05490 * std::sin(zm);
    sinzf = std::sin(zf);
    f2 = 0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * std::cos(zf);

    double sel = ee2_ * f2 + e3_ * f3;
    double sil = xi2_ * f2 + xi3_ * f3;
    double sll = xl2_ * f2 + xl3_ * f3 + xl4_ * sinzf;
    double sghl = xgh2_ * f2 + xgh3_ * f3 + xgh4_ * sinzf;
    double shll = xh2_ * f2 + xh3_ * f3;

    double pe = ses + sel;
    double pinc = sis + sil;
    double pl = sls + sll;
    double pgh = sghs + sghl;
    double ph = shs + shll;

    if (std::abs(inclo_) >= 0.2) {
        ph /= std::sin(inclo_);
        inclm += pinc;
        nodem += ph;
        argpm -= pgh;
    } else {
        // Apply Lyddane modification
        double siniq = std::sin(inclo_);
        double cosiq = std::cos(inclo_);
        double temp_mod = ph * cosiq;

        inclm += pinc;
        nodem += ph / siniq;
        argpm -= temp_mod / siniq;
    }

    em += pe;
    mm += pl;

    // Handle resonance effects
    if (irez_ != 0) {
        double ft = 0.0;
        constexpr double earthRotRate = 7.29211514668855e-5;  // Earth rotation rate rad/min

        // Synchronous resonance terms
        if (irez_ == 1) {
            double xndt;
            double xfact = mdot_ + argpdot_ + nodedot_ - earthRotRate - no_unkozai_;
            double xldot = xni_ + xfact;

            // Calculate resonance effects
            double theta = std::fmod(gsto_ + t * earthRotRate, SGP4_TWO_PI);
            xndt = del1_ * std::sin(xlamo_ - 2.0 * (nodeo_ + argpo_) + theta) +
                   del2_ * std::sin(2.0 * (xlamo_ - nodeo_ - argpo_)) +
                   del3_ * std::sin(3.0 * xlamo_ - nodeo_ - argpo_ + theta);

            if (std::abs(t - atime_) >= step) {
                // Integrate using predictor-corrector
                double stepp = step;
                if (t < atime_) {
                    stepp = -step;
                }

                while (std::abs(t - atime_) >= step) {
                    xldot = xni_ + xfact;
                    xli_ += xldot * stepp + xndt * step2;
                    xni_ += xndt * stepp;
                    atime_ += stepp;

                    theta = std::fmod(gsto_ + atime_ * earthRotRate, SGP4_TWO_PI);
                    xndt = del1_ * std::sin(xli_ - 2.0 * (nodeo_ + argpo_) + theta) +
                           del2_ * std::sin(2.0 * (xli_ - nodeo_ - argpo_)) +
                           del3_ * std::sin(3.0 * xli_ - nodeo_ - argpo_ + theta);
                }
            }

            ft = t - atime_;
            xldot = xni_ + xfact;
            nm = xni_ + xndt * ft;
            double xl = xli_ + xldot * ft + xndt * ft * ft * 0.5;
            mm = xl - 2.0 * nodem + 2.0 * theta;
        }

        // Half-day resonance terms
        if (irez_ == 2) {
            double xndt, xnddt, xomi;
            double xfact = mdot_ + dmdt_ + 2.0 * (nodedot_ + dnodt_ - earthRotRate) - no_unkozai_;
            double xldot = xni_ + xfact;
            constexpr double g22 = 5.7686396;
            constexpr double g32 = 0.95240898;
            constexpr double g44 = 1.8014998;
            constexpr double g52 = 1.0508330;
            constexpr double g54 = 4.4108898;

            double theta = std::fmod(gsto_ + t * earthRotRate, SGP4_TWO_PI);

            xomi = argpo_ + argpdot_ * atime_;
            double x2omi = xomi + xomi;
            double x2li = xli_ + xli_;

            xndt = d2201_ * std::sin(x2omi + xli_ - g22) +
                   d2211_ * std::sin(xli_ - g22) +
                   d3210_ * std::sin(xomi + xli_ - g32) +
                   d3222_ * std::sin(-xomi + xli_ - g32) +
                   d4410_ * std::sin(x2omi + x2li - g44) +
                   d4422_ * std::sin(x2li - g44) +
                   d5220_ * std::sin(xomi + xli_ - g52) +
                   d5232_ * std::sin(-xomi + xli_ - g52) +
                   d5421_ * std::sin(xomi + x2li - g54) +
                   d5433_ * std::sin(-xomi + x2li - g54);

            xnddt = d2201_ * std::cos(x2omi + xli_ - g22) +
                    d2211_ * std::cos(xli_ - g22) +
                    d3210_ * std::cos(xomi + xli_ - g32) +
                    d3222_ * std::cos(-xomi + xli_ - g32) +
                    d5220_ * std::cos(xomi + xli_ - g52) +
                    d5232_ * std::cos(-xomi + xli_ - g52) +
                    2.0 * (d4410_ * std::cos(x2omi + x2li - g44) +
                           d4422_ * std::cos(x2li - g44) +
                           d5421_ * std::cos(xomi + x2li - g54) +
                           d5433_ * std::cos(-xomi + x2li - g54));
            xnddt *= xldot;

            if (std::abs(t - atime_) >= step) {
                double stepp = step;
                if (t < atime_) {
                    stepp = -step;
                }

                while (std::abs(t - atime_) >= step) {
                    xldot = xni_ + xfact;
                    xli_ += xldot * stepp + xndt * step2;
                    xni_ += xndt * stepp;
                    atime_ += stepp;

                    xomi = argpo_ + argpdot_ * atime_;
                    x2omi = xomi + xomi;
                    x2li = xli_ + xli_;

                    xndt = d2201_ * std::sin(x2omi + xli_ - g22) +
                           d2211_ * std::sin(xli_ - g22) +
                           d3210_ * std::sin(xomi + xli_ - g32) +
                           d3222_ * std::sin(-xomi + xli_ - g32) +
                           d4410_ * std::sin(x2omi + x2li - g44) +
                           d4422_ * std::sin(x2li - g44) +
                           d5220_ * std::sin(xomi + xli_ - g52) +
                           d5232_ * std::sin(-xomi + xli_ - g52) +
                           d5421_ * std::sin(xomi + x2li - g54) +
                           d5433_ * std::sin(-xomi + x2li - g54);

                    xnddt = d2201_ * std::cos(x2omi + xli_ - g22) +
                            d2211_ * std::cos(xli_ - g22) +
                            d3210_ * std::cos(xomi + xli_ - g32) +
                            d3222_ * std::cos(-xomi + xli_ - g32) +
                            d5220_ * std::cos(xomi + xli_ - g52) +
                            d5232_ * std::cos(-xomi + xli_ - g52) +
                            2.0 * (d4410_ * std::cos(x2omi + x2li - g44) +
                                   d4422_ * std::cos(x2li - g44) +
                                   d5421_ * std::cos(xomi + x2li - g54) +
                                   d5433_ * std::cos(-xomi + x2li - g54));
                    xnddt *= xldot;
                }
            }

            ft = t - atime_;
            xldot = xni_ + xfact;
            nm = xni_ + xndt * ft + xnddt * ft * ft * 0.5;
            double xl = xli_ + xldot * ft + xndt * ft * ft * 0.5;
            double temp_mm = -nodem - nodem + theta + theta;
            mm = xl - xomi + temp_mm;
        }
    }
}

// Apply deep space periodic effects
void Orbit::deepSpacePeriodic(
    double t, double& em, double& inclm, double& nodem,
    double& argpm, double& mm) const {

    // Calculate lunar-solar periodics
    double zm = zmos_ + 0.017201977 * t;
    double zf = zm + 2.0 * 0.01675 * std::sin(zm);
    double sinzf = std::sin(zf);
    double f2 = 0.5 * sinzf * sinzf - 0.25;
    double f3 = -0.5 * sinzf * std::cos(zf);

    double ses = se2_ * f2 + se3_ * f3;
    double sis = si2_ * f2 + si3_ * f3;
    double sls = sl2_ * f2 + sl3_ * f3 + sl4_ * sinzf;
    double sghs = sgh2_ * f2 + sgh3_ * f3 + sgh4_ * sinzf;
    double shs = sh2_ * f2 + sh3_ * f3;

    zm = zmol_ + 0.22997150 * t;
    zf = zm + 2.0 * 0.05490 * std::sin(zm);
    sinzf = std::sin(zf);
    f2 = 0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * std::cos(zf);

    double sel = ee2_ * f2 + e3_ * f3;
    double sil = xi2_ * f2 + xi3_ * f3;
    double sll = xl2_ * f2 + xl3_ * f3 + xl4_ * sinzf;
    double sghl = xgh2_ * f2 + xgh3_ * f3 + xgh4_ * sinzf;
    double shll = xh2_ * f2 + xh3_ * f3;

    double pe = ses + sel - peo_;
    double pinc = sis + sil - pinco_;
    double pl = sls + sll - plo_;
    double pgh = sghs + sghl - pgho_;
    double ph = shs + shll - pho_;

    if (std::abs(inclo_) >= 0.2) {
        ph /= std::sin(inclo_);
        inclm += pinc;
        em += pe;
        nodem += ph;
        argpm -= pgh;
        mm += pl;
    } else {
        // Apply Lyddane modification
        double siniq = std::sin(inclo_);
        double cosiq = std::cos(inclo_);

        inclm += pinc;
        em += pe;

        double sinis = std::sin(inclm);
        double cosis = std::cos(inclm);

        if (std::abs(inclm) >= 0.2) {
            double temp_per = ph / sinis;
            nodem += temp_per;
            argpm -= pgh - cosiq * temp_per;
            mm += pl;
        } else {
            // Near-equatorial orbit, Lyddane modification
            double temp_per = ph * cosiq;
            nodem += ph / siniq;
            argpm -= temp_per / siniq;
            mm += pl;
        }
    }
}

// Core SGP4/SDP4 propagation
void Orbit::propagateSGP4(double tsince, double r[3], double v[3]) const {
    constexpr double radiusearthkm = SGP4_RADIUS_EARTH_KM;
    constexpr double xke = SGP4_XKE;
    constexpr double j2 = SGP4_J2;
    constexpr double j3oj2 = SGP4_J3OJ2;
    constexpr double vkmpersec = SGP4_VKMPERSEC;

    double cosio = std::cos(inclo_);
    double sinio = std::sin(inclo_);

    // Update for secular gravity and atmospheric drag
    double xmdf = mo_ + mdot_ * tsince;
    double argpdf = argpo_ + argpdot_ * tsince;
    double nodedf = nodeo_ + nodedot_ * tsince;
    double argpm = argpdf;
    double mm = xmdf;
    double t2 = tsince * tsince;
    double nodem = nodedf + nodecf_ * t2;
    double tempa = 1.0 - cc1_ * tsince;
    double tempe = bstar_ * cc4_ * tsince;
    double templ = t2cof_ * t2;

    if (!isimp_) {
        double delomg = omgcof_ * tsince;
        double delm = xmcof_ * (std::pow(1.0 + eta_ * std::cos(xmdf), 3) - delmo_);
        double temp_sgp = delomg + delm;
        mm = xmdf + temp_sgp;
        argpm = argpdf - temp_sgp;
        double t3 = t2 * tsince;
        double t4 = t3 * tsince;
        tempa = tempa - d2_ * t2 - d3_ * t3 - d4_ * t4;
        tempe = tempe + bstar_ * cc5_ * (std::sin(mm) - sinmao_);
        templ = templ + t3cof_ * t3 + t4 * (t4cof_ + tsince * t5cof_);
    }

    double nm = no_unkozai_;
    double em = ecco_;
    double inclm = inclo_;

    // Handle deep space satellites
    if (method_ == 'd') {
        deepSpaceSecular(tsince, em, argpm, inclm, nodem, mm, nm);
    }

    double am = std::pow(xke / nm, SGP4_X2O3) * tempa * tempa;
    nm = xke / std::pow(am, 1.5);
    em = em - tempe;

    // Check for eccentricity out of range
    if (em >= 1.0 || em < -0.001) {
        throw InvalidOrbitException("Eccentricity out of range during propagation: " + std::to_string(em));
    }
    if (em < 1.0e-6) {
        em = 1.0e-6;
    }

    mm = mm + no_unkozai_ * templ;
    double xlm = mm + argpm + nodem;

    nodem = std::fmod(nodem, SGP4_TWO_PI);
    argpm = std::fmod(argpm, SGP4_TWO_PI);
    xlm = std::fmod(xlm, SGP4_TWO_PI);
    mm = std::fmod(xlm - argpm - nodem, SGP4_TWO_PI);

    // Apply deep space periodic effects
    if (method_ == 'd') {
        deepSpacePeriodic(tsince, em, inclm, nodem, argpm, mm);
    }

    // Re-compute sini/cosi if inclination changed
    if (inclm != inclo_) {
        sinio = std::sin(inclm);
        cosio = std::cos(inclm);
    }

    // Check for eccentricity out of range (after deep space effects)
    if (em < 0.0) {
        em = 1.0e-6;
    }

    double sinim = std::sin(inclm);
    double cosim = std::cos(inclm);

    // Add lunar-solar periodics
    double ep = em;
    double xincp = inclm;
    double argpp = argpm;
    double nodep = nodem;
    double mp = mm;

    // Check if satellite decayed
    if (nm <= 0.0) {
        throw SatelliteDecayedException();
    }

    double eccsq = ep * ep;
    double omeosq = 1.0 - eccsq;

    if (omeosq <= 0.0) {
        throw InvalidOrbitException("Semi-latus rectum is negative");
    }

    double rteosq = std::sqrt(omeosq);

    // Long period periodics
    cosio = cosim;
    sinio = sinim;
    double cosio2 = cosio * cosio;

    double axnl = ep * std::cos(argpp);
    double temp_lp = 1.0 / (am * omeosq);
    double aynl = ep * std::sin(argpp) + temp_lp * aycof_;
    double xl = mp + argpp + nodep + temp_lp * xlcof_ * axnl;

    // Solve Kepler's equation
    double u = std::fmod(xl - nodep, SGP4_TWO_PI);
    double eo1 = u;
    double tem5 = 9999.9;
    int ktr = 1;
    double sineo1, coseo1;

    while ((std::abs(tem5) >= 1.0e-12) && (ktr <= 10)) {
        sineo1 = std::sin(eo1);
        coseo1 = std::cos(eo1);
        tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
        tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        if (std::abs(tem5) >= 0.95) {
            tem5 = tem5 > 0.0 ? 0.95 : -0.95;
        }
        eo1 = eo1 + tem5;
        ktr++;
    }

    // Short period preliminary quantities
    double ecose = axnl * coseo1 + aynl * sineo1;
    double esine = axnl * sineo1 - aynl * coseo1;
    double el2 = axnl * axnl + aynl * aynl;
    double pl = am * (1.0 - el2);

    if (pl < 0.0) {
        throw InvalidOrbitException("Semi-latus rectum is negative");
    }

    double rl = am * (1.0 - ecose);
    double rdotl = std::sqrt(am) * esine / rl;
    double rvdotl = std::sqrt(pl) / rl;
    double betal = std::sqrt(1.0 - el2);
    double temp_sp = esine / (1.0 + betal);
    double sinu = am / rl * (sineo1 - aynl - axnl * temp_sp);
    double cosu = am / rl * (coseo1 - axnl + aynl * temp_sp);
    double su = std::atan2(sinu, cosu);
    double sin2u = (cosu + cosu) * sinu;
    double cos2u = 1.0 - 2.0 * sinu * sinu;
    double temp_sp2 = 1.0 / pl;
    double temp1_sp = 0.5 * j2 * temp_sp2;
    double temp2_sp = temp1_sp * temp_sp2;

    // Update for short period periodics
    double con41 = 3.0 * cosio2 - 1.0;
    double x1mth2 = 1.0 - cosio2;
    double x7thm1 = 7.0 * cosio2 - 1.0;

    double mrt = rl * (1.0 - 1.5 * temp2_sp * betal * con41) + 0.5 * temp1_sp * x1mth2 * cos2u;
    su = su - 0.25 * temp2_sp * x7thm1 * sin2u;
    double xnode = nodep + 1.5 * temp2_sp * cosio * sin2u;
    double xinc = xincp + 1.5 * temp2_sp * cosio * sinio * cos2u;
    double mvt = rdotl - nm * temp1_sp * x1mth2 * sin2u / xke;
    double rvdot = rvdotl + nm * temp1_sp * (x1mth2 * cos2u + 1.5 * con41) / xke;

    // Orientation vectors
    double sinsu = std::sin(su);
    double cossu = std::cos(su);
    double snod = std::sin(xnode);
    double cnod = std::cos(xnode);
    double sini = std::sin(xinc);
    double cosi = std::cos(xinc);
    double xmx = -snod * cosi;
    double xmy = cnod * cosi;
    double ux = xmx * sinsu + cnod * cossu;
    double uy = xmy * sinsu + snod * cossu;
    double uz = sini * sinsu;
    double vx = xmx * cossu - cnod * sinsu;
    double vy = xmy * cossu - snod * sinsu;
    double vz = sini * cossu;

    // Position and velocity (in km and km/s)
    r[0] = mrt * ux * radiusearthkm;
    r[1] = mrt * uy * radiusearthkm;
    r[2] = mrt * uz * radiusearthkm;
    v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
    v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
    v[2] = (mvt * uz + rvdot * vz) * vkmpersec;

    // Check if satellite decayed
    if (mrt < 1.0) {
        throw SatelliteDecayedException();
    }
}

// Get position in ECI coordinates at a given Julian Date
Vec3 Orbit::getECI(double julianDate) const {
    if (!sgp4Initialized_) {
        initializeSGP4();
    }

    // Time since epoch in minutes
    double tsince = (julianDate - jdsatepoch_ - jdsatepochF_) * 1440.0;

    double r[3], v[3];
    propagateSGP4(tsince, r, v);

    return {r[0], r[1], r[2]};
}

// Get velocity in ECI coordinates at a given Julian Date
Vec3 Orbit::getVelocity(double julianDate) const {
    if (!sgp4Initialized_) {
        initializeSGP4();
    }

    // Time since epoch in minutes
    double tsince = (julianDate - jdsatepoch_ - jdsatepochF_) * 1440.0;

    double r[3], v[3];
    propagateSGP4(tsince, r, v);

    return {v[0], v[1], v[2]};
}

// Convert Earth Centered Inertial (ECI) to Earth Centered Earth Fixed (ECEF) using Greenwich Sidereal Time
Vec3 eciToECEF(const Vec3 &eci, double gst) {
    double cosGST = std::cos(gst);
    double sinGST = std::sin(gst);
    
    return {
         eci.x * cosGST + eci.y * sinGST,
        -eci.x * sinGST + eci.y * cosGST,
         eci.z
    };
}

// WGS84 ellipsoid constants (used by multiple functions)
constexpr double WGS84_A = 6378.137;              // Semi-major axis (km) - equatorial radius
constexpr double WGS84_F = 1.0 / 298.257223563;   // Flattening
constexpr double WGS84_E2 = WGS84_F * (2 - WGS84_F);  // Eccentricity squared ≈ 0.00669437999014

// Convert ECEF coordinates to geodetic latitude, longitude, and altitude
Geodetic ecefToGeodetic(const Vec3 &ecef) {
    double x = ecef.x, y = ecef.y, z = ecef.z;
    double lon = std::atan2(y, x);
    double p = std::sqrt(x*x + y*y);

    // Iterative latitude calculation (Bowring's method)
    double lat = std::atan2(z, p * (1 - WGS84_E2));  // initial guess
    for (int i = 0; i < 10; ++i) {
        double sinLat = std::sin(lat);
        double N = WGS84_A / std::sqrt(1 - WGS84_E2 * sinLat * sinLat);
        lat = std::atan2(z + WGS84_E2 * N * sinLat, p);
    }

    double sinLat = std::sin(lat);
    double N = WGS84_A / std::sqrt(1 - WGS84_E2 * sinLat * sinLat);
    double alt = p / std::cos(lat) - N;

    return {lat, lon, alt};
}

// Convert geodetic coordinates to ECEF (Earth-Centered Earth-Fixed)
// See documentation in orbit.hpp for detailed algorithm explanation
Vec3 Geodetic::toECEF() const {
    double sinLat = std::sin(latInRadians);
    double cosLat = std::cos(latInRadians);
    double sinLon = std::sin(lonInRadians);
    double cosLon = std::cos(lonInRadians);

    // Radius of curvature in the prime vertical
    // This is the distance from the surface to the Z-axis along the ellipsoid normal
    double N = WGS84_A / std::sqrt(1.0 - WGS84_E2 * sinLat * sinLat);

    // ECEF coordinates
    // The (1 - e²) factor in Z accounts for the ellipsoid's polar flattening
    return {
        (N + altInKilometers) * cosLat * cosLon,
        (N + altInKilometers) * cosLat * sinLon,
        (N * (1.0 - WGS84_E2) + altInKilometers) * sinLat
    };
}

// Transform ECEF coordinates to ENU (East-North-Up) local tangent plane
// See documentation in orbit.hpp for detailed algorithm explanation
Vec3 ecefToENU(const Vec3& targetECEF, const Geodetic& observer) {
    // Step 1: Get observer's ECEF position and compute difference vector
    Vec3 observerECEF = observer.toECEF();
    Vec3 diff = targetECEF - observerECEF;

    // Step 2: Precompute trig values for rotation matrix
    double sinLat = std::sin(observer.latInRadians);
    double cosLat = std::cos(observer.latInRadians);
    double sinLon = std::sin(observer.lonInRadians);
    double cosLon = std::cos(observer.lonInRadians);

    // Step 3: Apply rotation matrix (ECEF to ENU)
    // This matrix rotates the ECEF frame to align with the local horizon:
    //   - East points along the local latitude circle (toward increasing longitude)
    //   - North points along the local meridian (toward the pole)
    //   - Up points radially outward (normal to the ellipsoid)
    double east  = -sinLon * diff.x + cosLon * diff.y;
    double north = -sinLat * cosLon * diff.x - sinLat * sinLon * diff.y + cosLat * diff.z;
    double up    =  cosLat * cosLon * diff.x + cosLat * sinLon * diff.y + sinLat * diff.z;

    return {east, north, up};
}

// Compute look angles from observer to a target in ECEF coordinates
LookAngles getLookAngles(const Vec3& satECEF, const Geodetic& observer) {
    // Transform to local ENU coordinates
    Vec3 enu = ecefToENU(satECEF, observer);

    // Compute slant range (straight-line distance)
    double range = enu.magnitude();

    // Compute elevation angle
    // elevation = arcsin(Up / range)
    // When satellite is directly overhead, Up = range, so elevation = π/2
    // When satellite is on horizon, Up = 0, so elevation = 0
    double elevation = std::asin(enu.z / range);

    // Compute azimuth angle
    // azimuth = arctan2(East, North)
    // This gives: 0 = North, π/2 = East, π = South, -π/2 = West
    double azimuth = std::atan2(enu.x, enu.y);

    // Normalize azimuth to [0, 2π)
    if (azimuth < 0) {
        azimuth += 2.0 * std::numbers::pi;
    }

    return {azimuth, elevation, range};
}

// Compute look angles from observer to a satellite at a specific time
LookAngles getLookAngles(const Orbit& orbit, const Geodetic& observer, double julianDate) {
    // Get ECI coordinates using SGP4
    Vec3 eci = orbit.getECI(julianDate);

    // Convert ECI to ECEF using Greenwich Sidereal Time
    double gst = gmst(julianDate);
    Vec3 ecef = eciToECEF(eci, gst);

    // Compute look angles
    return getLookAngles(ecef, observer);
}

// Check visibility based on look angles
bool isVisible(const LookAngles& angles, double minElevationInRadians) {
    return angles.elevationInRadians >= minElevationInRadians;
}

// Check visibility of a satellite at a specific time
bool isVisible(const Orbit& orbit, const Geodetic& observer,
               double julianDate, double minElevationInRadians) {
    LookAngles angles = getLookAngles(orbit, observer, julianDate);
    return isVisible(angles, minElevationInRadians);
}

// Helper: Get elevation at a specific time
static double getElevationAtTime(const Orbit& orbit, const Geodetic& observer, time_point t) {
    double jd = toJulianDate(t);
    LookAngles angles = getLookAngles(orbit, observer, jd);
    return angles.elevationInRadians;
}

// Helper: Binary search to find precise threshold crossing time
// Returns the time when elevation crosses the threshold
// 'rising' indicates if we're looking for a rise (true) or set (false)
static time_point binarySearchCrossing(
    const Orbit& orbit,
    const Geodetic& observer,
    double threshold,
    time_point low,
    time_point high,
    bool rising) {

    using namespace std::chrono;

    // Refine to within 1 second
    while (duration_cast<seconds>(high - low).count() > 1) {
        time_point mid = low + (high - low) / 2;
        double elev = getElevationAtTime(orbit, observer, mid);

        bool aboveThreshold = elev >= threshold;

        if (rising) {
            // Looking for rise: below threshold -> above threshold
            // If mid is above, the crossing is before mid
            if (aboveThreshold) {
                high = mid;
            } else {
                low = mid;
            }
        } else {
            // Looking for set: above threshold -> below threshold
            // If mid is above, the crossing is after mid
            if (aboveThreshold) {
                low = mid;
            } else {
                high = mid;
            }
        }
    }

    return rising ? high : low;
}

// Helper: Golden section search to find maximum elevation time
// Returns the time when elevation is at maximum between low and high
static time_point goldenSectionSearchMax(
    const Orbit& orbit,
    const Geodetic& observer,
    time_point low,
    time_point high) {

    using namespace std::chrono;

    constexpr double PHI = 1.618033988749895;  // Golden ratio
    constexpr double RESPHI = 2.0 - PHI;       // 1/phi

    auto duration = high - low;
    time_point x1 = low + duration_cast<system_clock::duration>(duration * RESPHI);
    time_point x2 = high - duration_cast<system_clock::duration>(duration * RESPHI);

    double f1 = getElevationAtTime(orbit, observer, x1);
    double f2 = getElevationAtTime(orbit, observer, x2);

    // Refine to within 1 second
    while (duration_cast<seconds>(high - low).count() > 1) {
        if (f1 > f2) {
            high = x2;
            x2 = x1;
            f2 = f1;
            duration = high - low;
            x1 = low + duration_cast<system_clock::duration>(duration * RESPHI);
            f1 = getElevationAtTime(orbit, observer, x1);
        } else {
            low = x1;
            x1 = x2;
            f1 = f2;
            duration = high - low;
            x2 = high - duration_cast<system_clock::duration>(duration * RESPHI);
            f2 = getElevationAtTime(orbit, observer, x2);
        }
    }

    // Return the midpoint
    return low + (high - low) / 2;
}

// Find the next satellite pass over an observer's location
std::optional<PassInfo> findNextPass(
    const Orbit& orbit,
    const Geodetic& observer,
    double minElevationInRadians,
    time_point startTime) {

    using namespace std::chrono;

    // Search parameters
    constexpr auto COARSE_STEP = seconds(60);        // 60-second steps for coarse search
    constexpr auto MAX_SEARCH = hours(48);           // Search up to 48 hours ahead

    auto searchEnd = startTime + MAX_SEARCH;
    auto currentTime = startTime;

    // Track visibility state
    bool wasVisible = isVisible(orbit, observer, toJulianDate(currentTime), minElevationInRadians);

    // If we start during a pass, skip to the end of it first
    if (wasVisible) {
        while (currentTime < searchEnd) {
            currentTime += COARSE_STEP;
            bool nowVisible = isVisible(orbit, observer, toJulianDate(currentTime), minElevationInRadians);
            if (!nowVisible) {
                wasVisible = false;
                break;
            }
        }
        if (wasVisible) {
            // Still visible after 48 hours - probably a GEO satellite or error
            return std::nullopt;
        }
    }

    // Coarse search: find when satellite rises above threshold
    time_point coarseRise, coarseSet;
    bool foundRise = false;
    bool foundSet = false;

    while (currentTime < searchEnd) {
        currentTime += COARSE_STEP;
        bool nowVisible = isVisible(orbit, observer, toJulianDate(currentTime), minElevationInRadians);

        if (!wasVisible && nowVisible) {
            // Found rise
            coarseRise = currentTime - COARSE_STEP;
            foundRise = true;
            wasVisible = true;
        } else if (wasVisible && !nowVisible) {
            // Found set
            coarseSet = currentTime;
            foundSet = true;
            break;
        }

        wasVisible = nowVisible;
    }

    // Handle edge cases
    if (!foundRise) {
        return std::nullopt;  // No pass found in search window
    }

    if (!foundSet) {
        // Pass started but didn't end - use search end as approximate set
        coarseSet = searchEnd;
    }

    // Binary search refinement for precise rise time
    time_point preciseRise = binarySearchCrossing(
        orbit, observer, minElevationInRadians,
        coarseRise, coarseRise + COARSE_STEP, true);

    // Binary search refinement for precise set time
    time_point preciseSet = binarySearchCrossing(
        orbit, observer, minElevationInRadians,
        coarseSet - COARSE_STEP, coarseSet, false);

    // Golden section search for maximum elevation time
    time_point maxTime = goldenSectionSearchMax(orbit, observer, preciseRise, preciseSet);

    // Compute look angles at each key time
    double riseJD = toJulianDate(preciseRise);
    double maxJD = toJulianDate(maxTime);
    double setJD = toJulianDate(preciseSet);

    LookAngles riseAngles = getLookAngles(orbit, observer, riseJD);
    LookAngles maxAngles = getLookAngles(orbit, observer, maxJD);
    LookAngles setAngles = getLookAngles(orbit, observer, setJD);

    return PassInfo{
        .noradID = orbit.getNoradID(),
        .name = orbit.getName(),
        .riseTime = preciseRise,
        .maxElevationTime = maxTime,
        .setTime = preciseSet,
        .riseAngles = riseAngles,
        .maxAngles = maxAngles,
        .setAngles = setAngles
    };
}

// Get geodetic location (lat, lon, alt) of the satellite at a given time
Geodetic Orbit::getGeodeticLocationAtTime(const time_point tp) const {
    double julianDate = toJulianDate(tp);
    return getGeodeticLocationAtTime(julianDate);
}   

// Get geodetic location (lat, lon, alt) of the satellite at a given time
Geodetic Orbit::getGeodeticLocationAtTime(const double julianDate) const {
    // Get ECI coordinates using SGP4
    Vec3 eci = getECI(julianDate);

    // Get GMST at the requested time (the earth's rotation angle)
    double gst = gmst(julianDate);

    // Convert ECI to ECEF (the satellite's position in earth-fixed frame)
    Vec3 ecef = eciToECEF(eci, gst);

    // Convert ECEF to geodetic coordinates (latitude/longitude/altitude)
    return ecefToGeodetic(ecef);
}

// Update orbital elements from TLE data with a seperate name
void Orbit::updateFromTLE(const std::string_view &name, const std::string_view &tle) {
    updateFromTLE(tle);
    this->name = std::string(name);
}

// Update orbital elements from TLE data
void Orbit::updateFromTLE(const std::string_view &tle) {
    // Reset SGP4 initialization state so it re-initializes on next propagation
    sgp4Initialized_ = false;

    bool firstLineParsed = false;
    bool secondLineParsed = false;
    for (auto line : tle | std::views::split('\n')) {
        std::string_view lineView(line.begin(), line.end());
        lineView = trimLeft(trimRight(lineView));
        if (lineView.starts_with("1 ")) {
            // NORAD ID is columns 3-7
            noradID = toNumber<int>(trimLeft(lineView.substr(2, 5)));
            // Classification is column 8
            classification = lineView[7];
            // Designator is columns 10-17
            designator = std::string(lineView.substr(9, 8));
            // First Derivative of Mean Motion is columns 34-43
            firstDerivativeMeanMotion = toNumber<double>(trimLeft(lineView.substr(33, 10)));
            // Second Derivative of Mean Motion is columns 45-52 (exponential format)
            secondDerivativeMeanMotion = fromExponentialString(trimLeft(lineView.substr(44, 8)));
            // Bstar Drag Term is columns 54-61 (exponential format)
            bstarDragTerm = fromExponentialString(trimLeft(lineView.substr(53, 8)));
            // Time since epoch is columns 19-32
            epoch = parseEpoch(lineView.substr(18, 14));
            // Element Set Number is columns 65-68
            elementSetNumber = toNumber<int>(trimLeft(lineView.substr(64, 4)));
            firstLineParsed = true;
        } else if (lineView.starts_with("2 ")) {
            // Inclination is columns 9-16
            inclination = toNumber<double>(trimLeft(lineView.substr(8, 8)));
            // RAAN is columns 18-25
            rightAscensionOfAscendingNode = toNumber<double>(trimLeft(lineView.substr(17, 8)));
            // Eccentricity is columns 27-33 (decimal implied)
            eccentricity = toNumber<double>("0." + std::string(trimLeft(lineView.substr(26, 7))));
            // Argument of perigee is columns 35-42
            argumentOfPerigee = toNumber<double>(trimLeft(lineView.substr(34, 8)));
            // Mean Anomaly is columns 44-51
            meanAnomaly = toNumber<double>(trimLeft(lineView.substr(43, 8)));
            // Mean Motion is columns 53-63
            meanMotion = toNumber<double>(trimLeft(lineView.substr(52, 11)));
            // Revolution number at epoch is columns 64-68
            revolutionNumberAtEpoch = toNumber<int>(trimLeft(lineView.substr(63, 5)));
            secondLineParsed = true;
        } else if (!firstLineParsed && !secondLineParsed && !lineView.empty()) {
            // If neither line has been parsed and this line is not empty, assume it's the name
            name = std::string(lineView);
        }

        if (firstLineParsed && secondLineParsed) {
            break;
        }
    }
}

int Orbit::getNoradID() const {
    return noradID;
}

char Orbit::getClassification() const {
    return classification;
}

std::string Orbit::getDesignator() const {
    return designator;
}

time_point Orbit::getEpoch() const {
    return epoch;
}

double Orbit::getFirstDerivativeMeanMotion() const {
    return firstDerivativeMeanMotion;
}

double Orbit::getSecondDerivativeMeanMotion() const {
    return secondDerivativeMeanMotion;
}

double Orbit::getBstarDragTerm() const {
    return bstarDragTerm;
}

int Orbit::getElementSetNumber() const {
    return elementSetNumber;
}

double Orbit::getInclination() const {
    return inclination;
}

double Orbit::getRightAscensionOfAscendingNode() const {
    return rightAscensionOfAscendingNode;
}

double Orbit::getEccentricity() const {
    return eccentricity;
}

double Orbit::getArgumentOfPerigee() const {
    return argumentOfPerigee;
}

double Orbit::getMeanAnomaly() const {
    return meanAnomaly;
}

double Orbit::getMeanMotion() const {
    return meanMotion;
}

int Orbit::getRevolutionNumberAtEpoch() const {
    return revolutionNumberAtEpoch;
}

std::string Orbit::getName() const {
    return name;
}

void Orbit::printInfo(std::ostream &os) const {
    os << getName() << std::endl;
    os << "  NORAD ID: " << getNoradID() << std::endl;
    os << "  Classification: " << getClassification() << std::endl;
    os << "  Designator: " << getDesignator() << std::endl;
    auto epoch = std::chrono::floor<std::chrono::seconds>(getEpoch());
    os << "  Epoch: " << date::format("%F %T UTC", epoch) << std::endl;
    os << "  First Derivative of Mean Motion: " << getFirstDerivativeMeanMotion() << std::endl;
    os << "  Second Derivative of Mean Motion: " << getSecondDerivativeMeanMotion() << std::endl;
    os << "  Bstar Drag Term: " << getBstarDragTerm() << std::endl;
    os << "  Element Set Number: " << getElementSetNumber() << std::endl;
    os << "  Inclination: " << getInclination() << " deg" << std::endl;
    os << "  Right Ascension of Ascending Node: " << getRightAscensionOfAscendingNode() << " deg" << std::endl;
    os << "  Eccentricity: " << getEccentricity() << std::endl;
    os << "  Argument of Perigee: " << getArgumentOfPerigee() << " deg" << std::endl;
    os << "  Mean Anomaly: " << getMeanAnomaly() << " deg" << std::endl;
    os << "  Mean Motion: " << getMeanMotion() << " revs per day" << std::endl;
    os << "  Revolution Number at Epoch: " << getRevolutionNumberAtEpoch() << std::endl;
    os << std::endl;
}

// Calculate TLE line checksum (mod 10 sum of digits, with '-' counting as 1)
int calculateChecksum(const std::string& line) {
    int sum = 0;
    for (char c : line) {
        if (c >= '0' && c <= '9') {
            sum += (c - '0');
        } else if (c == '-') {
            sum += 1;
        }
    }
    return sum % 10;
}

// Format a value in TLE exponential notation (e.g., " 00000+0" or " 15237-3" or "-12345-6")
// The format is: [sign]NNNNN[sign]E where NNNNN is 5 digits of mantissa and E is exponent
std::string toTLEExponential(double value) {
    if (value == 0.0) {
        return " 00000+0";
    }

    char sign = (value >= 0) ? ' ' : '-';
    value = std::abs(value);

    // Get exponent (power of 10)
    int exponent = static_cast<int>(std::floor(std::log10(value)));

    // Normalize mantissa to be in range [0.1, 1.0)
    double mantissa = value / std::pow(10.0, exponent + 1);

    // Convert mantissa to 5-digit integer
    int mantissaInt = static_cast<int>(std::round(mantissa * 100000));

    // Handle rounding overflow
    if (mantissaInt >= 100000) {
        mantissaInt = 10000;
        exponent++;
    }

    // Format exponent with sign
    char expSign = (exponent + 1 >= 0) ? '+' : '-';
    int expAbs = std::abs(exponent + 1);

    std::ostringstream ss;
    ss << sign << std::setw(5) << std::setfill('0') << mantissaInt << expSign << expAbs;
    return ss.str();
}

// Format first derivative of mean motion for TLE (e.g., " .00008010" or "-.00012345")
std::string formatFirstDerivative(double value) {
    char sign = (value >= 0) ? ' ' : '-';
    value = std::abs(value);

    std::ostringstream ss;
    ss << sign << '.' << std::setw(8) << std::setfill('0')
       << static_cast<int>(std::round(value * 100000000));
    return ss.str();
}

// Get standard 3-line TLE representation of the orbital elements
std::string Orbit::getTLE() const {
    using namespace std::chrono;

    std::ostringstream tleStream;

    // Line 0: Name
    tleStream << name << '\n';

    // Line 1
    // Compute epoch in TLE format: YYDDD.DDDDDDDD
    auto epochDays = floor<days>(epoch);
    year_month_day ymd{epochDays};
    int year = static_cast<int>(ymd.year());
    int twoDigitYear = year % 100;

    // Day of year (1-366)
    auto yearStart = sys_days{ymd.year()/January/1};
    int dayOfYear = (epochDays - yearStart).count() + 1;

    // Fractional part of day
    auto timeOfDay = epoch - epochDays;
    double fracDay = duration_cast<duration<double, std::ratio<86400>>>(timeOfDay).count();

    // Format epoch string (14 chars: YYDDD.DDDDDDDD)
    std::ostringstream epochStr;
    epochStr << std::setw(2) << std::setfill('0') << twoDigitYear
             << std::setw(3) << std::setfill('0') << dayOfYear
             << '.' << std::setw(8) << std::setfill('0')
             << static_cast<long>(std::round(fracDay * 100000000));

    // Build line 1 without checksum (68 chars), then add checksum
    // TLE Line 1 format (columns 1-indexed):
    // Col 01: Line number (1)
    // Col 02: Space
    // Col 03-07: NORAD Catalog Number (5 digits)
    // Col 08: Classification (U/C/S)
    // Col 09: Space
    // Col 10-17: International Designator (8 chars, left-justified, space-padded)
    // Col 18: Space
    // Col 19-32: Epoch (14 chars: YYDDD.DDDDDDDD)
    // Col 33: Space
    // Col 34-43: First Derivative Mean Motion (10 chars: sX.XXXXXXXX)
    // Col 44: Space
    // Col 45-52: Second Derivative Mean Motion (8 chars: sNNNNNsE)
    // Col 53: Space
    // Col 54-61: BSTAR (8 chars: sNNNNNsE)
    // Col 62: Space
    // Col 63: Ephemeris Type (always 0)
    // Col 64: Space
    // Col 65-68: Element Set Number (4 digits, right-justified, space-padded)
    // Col 69: Checksum
    std::ostringstream line1;
    line1 << "1 "
          << std::setw(5) << std::setfill('0') << noradID
          << classification << ' '
          << std::left << std::setw(8) << std::setfill(' ') << designator << ' '
          << epochStr.str() << ' '
          << formatFirstDerivative(firstDerivativeMeanMotion) << ' '
          << toTLEExponential(secondDerivativeMeanMotion) << ' '
          << toTLEExponential(bstarDragTerm) << ' '
          << "0 "
          << std::right << std::setw(4) << std::setfill(' ') << (elementSetNumber % 10000);

    std::string line1Str = line1.str();
    int checksum1 = calculateChecksum(line1Str);
    tleStream << line1Str << checksum1 << '\n';

    // Build line 2 without checksum (68 chars), then add checksum
    // TLE Line 2 format (columns 1-indexed):
    // Col 01: Line number (2)
    // Col 02: Space
    // Col 03-07: NORAD Catalog Number (5 digits)
    // Col 08: Space
    // Col 09-16: Inclination (8 chars: XXX.XXXX)
    // Col 17: Space
    // Col 18-25: RAAN (8 chars: XXX.XXXX)
    // Col 26: Space
    // Col 27-33: Eccentricity (7 digits, implied decimal)
    // Col 34: Space
    // Col 35-42: Argument of Perigee (8 chars: XXX.XXXX)
    // Col 43: Space
    // Col 44-51: Mean Anomaly (8 chars: XXX.XXXX)
    // Col 52: Space
    // Col 53-63: Mean Motion (11 chars: XX.XXXXXXXX)
    // Col 64-68: Revolution Number at Epoch (5 digits)
    // Col 69: Checksum
    int eccInt = static_cast<int>(std::round(eccentricity * 10000000));

    std::ostringstream line2;
    line2 << "2 "
          << std::setw(5) << std::setfill('0') << noradID << ' '
          << std::right << std::setw(8) << std::setfill(' ')
          << std::fixed << std::setprecision(4) << inclination << ' '
          << std::right << std::setw(8) << std::setfill(' ')
          << std::fixed << std::setprecision(4) << rightAscensionOfAscendingNode << ' '
          << std::setw(7) << std::setfill('0') << eccInt << ' '
          << std::right << std::setw(8) << std::setfill(' ')
          << std::fixed << std::setprecision(4) << argumentOfPerigee << ' '
          << std::right << std::setw(8) << std::setfill(' ')
          << std::fixed << std::setprecision(4) << meanAnomaly << ' '
          << std::right << std::setw(11) << std::setfill(' ')
          << std::fixed << std::setprecision(8) << meanMotion
          << std::setw(5) << std::setfill('0') << (revolutionNumberAtEpoch % 100000);

    std::string line2Str = line2.str();
    int checksum2 = calculateChecksum(line2Str);
    tleStream << line2Str << checksum2 << '\n';

    return tleStream.str();
}

// Load TLE database from file using the standard 3-line TLE format
void loadTLEDatabase(const std::string &filepath, std::map<int, Orbit> &database, std::ostream &logStream) {
    logStream << "Loading TLE database from file: " << filepath << std::endl;
    if (!std::filesystem::exists(filepath)) {
        logStream << "TLE database file does not exist: " + filepath << std::endl;
        return;
    }
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open TLE database file: " + filepath);
    }
    loadTLEDatabase(file, database, logStream);
}

// Load TLE database from stream using the standard 3-line TLE format
void loadTLEDatabase(std::istream &s, std::map<int, Orbit> &database, std::ostream &logStream) {
    bool haveFirstLine = false;
    bool haveSecondLine = false;
    std::string line, line1, line2, nameLine;
    int entriesLoaded = 0;
    logStream << "Loading TLE database... " << std::flush;
    while (std::getline(s, line)) {
        if (line.empty()) continue;

        if (line.starts_with("1 ")) {
            line1 = line;
            haveFirstLine = true;
        } else if (line.starts_with("2 ")) {
            line2 = line;
            haveSecondLine = true;
        } else {
            nameLine = line; // Assume it's the name
        }

        if (haveFirstLine && haveSecondLine) {
            Orbit orbit;

            std::ostringstream tleStream;
            tleStream << nameLine << '\n' << line1 << '\n' << line2 << '\n';

            orbit.updateFromTLE(tleStream.str());

            database[orbit.getNoradID()] = orbit;

            haveFirstLine = false;
            haveSecondLine = false;
            line1.clear();
            line2.clear();
            nameLine.clear();
            entriesLoaded++;
        }
    }
    logStream << " done." << std::endl;
    logStream << "Loaded " << entriesLoaded << " TLE entries." << std::endl;
    logStream << std::flush;
}

// Save TLE database to file using the standard 3-line TLE format
void saveTLEDatabase(const std::string &filepath, const std::map<int, Orbit> &database, std::ostream &logStream) {
    std::ofstream file(filepath);
    logStream << "Saving TLE database to file: " << filepath << std::endl;
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + filepath);
    }
    saveTLEDatabase(file, database, logStream);
}

// Save TLE database to stream using the standard 3-line TLE format
void saveTLEDatabase(std::ostream &s, const std::map<int, Orbit> &database, std::ostream &logStream) {
    logStream << "Saving TLE database... " << std::flush;
    int count = 0;
    for (const auto &[id, orbit] : database) {
        s << orbit.getTLE();
        count++;
    }
    logStream << " done." << std::endl;
    logStream << "Saved " << count << " TLE entries." << std::endl;
    logStream << std::flush;
}

} // namespace sattrack
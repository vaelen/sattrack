# SatTrack

SatTrack is a command-line utility for tracking satellites and predicting when they will be visible from a given location. It's targeted mainly at the amateur radio and amateur satellite communities.

It has the following features:

1. Download TLE data from the Internet
2. Display orbital elements for one or more satellites
3. Calculate the position of a satellite at a given time
4. Calculate look angles (azimuth/elevation/range) for antenna pointing
5. Check if a satellite is currently visible from your location
6. Real-time satellite tracking with configurable update interval
7. Predict upcoming satellite passes

## Build instructions

Install dependencies:

```sh
sudo apt install build-essential cmake libcurl4-openssl-dev libcurlpp-dev libspdlog-dev libfmt-dev libgtest-dev
```

Build the application:

```sh
mkdir build
cd build
cmake ..
make
make test
```

### Cross-compilation for ARM (BeagleBone, Raspberry Pi, etc.)

For a fully static build that runs on embedded ARM systems with older libraries:

**Install the musl cross-compiler:**

```sh
cd /opt
sudo wget https://musl.cc/arm-linux-musleabihf-cross.tgz
sudo tar xzf arm-linux-musleabihf-cross.tgz
sudo rm arm-linux-musleabihf-cross.tgz
```

Add the compiler to your PATH (add to `~/.bashrc` for persistence):

```sh
export PATH="/opt/arm-linux-musleabihf-cross/bin:$PATH"
```

**Build the ARM sysroot (one-time):**

This builds static versions of zlib, OpenSSL, curl, and curlpp for ARM:

```sh
TOOLCHAIN_PREFIX=arm-linux-musleabihf SYSROOT=~/arm-sysroot-musl ./scripts/build-arm-sysroot.sh
```

**Build for ARM:**

```sh
mkdir build-arm && cd build-arm
cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/arm-linux-musleabihf.cmake -DCMAKE_BUILD_TYPE=Release
make
```

**Strip the binary (optional, reduces size from ~17MB to ~5MB):**

```sh
arm-linux-musleabihf-strip sattrack
```

The resulting binary is fully static with no runtime dependencies and will run on any armv7l Linux system regardless of the installed libc version.

## Implementation details

SatTrack is a command-line application built using C++20. I wrote it largely as an opportunity to better familiarize myself with the features of C++11, C++14, C++17, and C++20 since I originally learned C++ back in the late 90s.

I used Claude Code to help speed up my implementation of SatTrack by helping me maintain my build scripts, write test cases, update documentation, and perform research on calculations and C++20 features that I wasn't familiar with, but most of the actual code was written by me.

Although I wrote the initial [Keplerian](https://en.wikipedia.org/wiki/Kepler_orbit) propagation code by hand and learned a lot from doing so, when I decided to move to [SGP4/SDP4](https://en.wikipedia.org/wiki/Simplified_perturbations_models) to improve pass prediction accuracy I had Claude port the SGP4/SDP4 [reference implementation](https://celestrak.org/software/vallado-sw.php) published by Celestrak rather than writing it myself. I didn't use an off-the-shelf library because I wanted tighter integration with my other code and with the C++20 standard library.

## Data sources

SatTrack uses [CelesTrak](https://celestrak.org/) as its primary source for TLE data. CelesTrak is a non-profit service that provides satellite information for free. If you feel so inclined, they take donations to help cover their costs.

## Configuration

All configuration options can be provided on the command line, but you can also create a config file to store common options, like the list of satellites you care about. The config file is stored in TOML format and the options available are the same as the command line options. The config file is located at `~/.sattrack.toml` by default, but you can also supply a path to the config file using the `--config` command line option.

Here is an example config file:

```toml
# Station Location
lat = 35.58
long = 139.48
alt = 30

info.id = 25544
tle.id = 25544
```

## Usage

The most common usage pattern looks like this:

```sh
# First update our TLE in case it is out of date
$ sattrack update

Loading TLE database from file: /home/andrew/.sattrack.tle
Loading TLE database...  done.
Loaded 13526 TLE entries.
Downloading TLE data for group: active... Done.
Saving TLE database to file: /home/andrew/.sattrack.tle
Saving TLE database...  done.
Saved 13526 TLE entries.

# Show upcomming passes
$ sattrack passes 25544

Passes for ISS (ZARYA):

Upcoming Passes:
        Satellite                   Start                      End               Duration        Starts In        Start Az      End Az       Max Az     Max Elev 
------------------------- ------------------------- ------------------------- -------------- ------------------ ------------ ------------ ------------ ----------
 25544 ISS (ZARYA)         2025-12-03 00:02:09 UTC   2025-12-03 00:06:01 UTC       3m 52s        0h 17m  9s      162.51 SSE    92.23 E     127.40 SE     14.56   
 25544 ISS (ZARYA)         2025-12-03 01:37:18 UTC   2025-12-03 01:43:52 UTC       6m 33s        1h 52m 18s      245.78 WSW    37.47 NE    322.04 NW     47.48   
 25544 ISS (ZARYA)         2025-12-03 03:17:46 UTC   2025-12-03 03:18:20 UTC       0m 33s        3h 32m 46s      334.48 NNW   343.97 NNW   339.22 NNW    10.08   
 25544 ISS (ZARYA)         2025-12-03 08:08:50 UTC   2025-12-03 08:15:20 UTC       6m 30s        8h 23m 50s      323.54 NW    112.49 ESE    37.93 NE     44.74   
 25544 ISS (ZARYA)         2025-12-03 09:46:33 UTC   2025-12-03 09:50:39 UTC       4m  5s       10h  1m 33s      269.94 W     194.63 SSW   232.24 SW     15.25   
 25544 ISS (ZARYA)         2025-12-04 00:49:48 UTC   2025-12-04 00:56:25 UTC       6m 37s       25h  4m 48s      213.26 SSW    56.44 ENE   135.35 SE     55.23  
```

### Commands

SatTrack provides the following commands:

#### `update` - Update local TLE data from Celestrak

Downloads TLE data from Celestrak and stores it locally. You can specify which satellite groups to download.

```sh
# Update with the default "active" group
sattrack update

# Update specific groups
sattrack update amateur visual
```

#### `info` - View satellite information

Displays orbital elements for satellites from the local TLE database.

```sh
sattrack info 25544
```

#### `tle` - View raw TLE data

Displays the raw Two-Line Element set data from the local TLE database.

```sh
sattrack tle 25544
```

#### `geo` - Get geodetic location

Calculates the geodetic location (latitude, longitude, altitude) of a satellite at a given time.

```sh
sattrack geo 25544 --time "2025-12-01 12:00:00"
```

#### `look` - Get look angles

Calculates the look angles (azimuth, elevation, range) for antenna pointing. Uses TLE data from the local database.

Options:

- `--time <time>` - Time at which to get look angles (format: YYYY-MM-DD HH:MM:SS UTC, default: now)

```sh
sattrack --lat 35.58 --long 139.48 look 25544
sattrack --lat 35.58 --long 139.48 look 25544 --time "2025-12-01 12:00:00"
```

#### `visible` - Check satellite visibility

Checks if a satellite is currently visible (above the minimum elevation threshold) from your location.

Options:

- `--elev <n>` - Minimum elevation above the horizon in degrees (default 10)
- `--time <time>` - Time at which to check visibility (format: YYYY-MM-DD HH:MM:SS UTC, default: now)

```sh
sattrack --lat 35.58 --long 139.48 visible 25544
sattrack --lat 35.58 --long 139.48 visible 25544 --elev 5
```

#### `track` - Real-time tracking

Provides real-time tracking output with configurable update interval and duration.

Options:

- `--interval <n>` - Update interval in seconds (default 1)
- `--duration <n>` - Duration to track in seconds (default 60)

```sh
sattrack --lat 35.58 --long 139.48 track 25544
sattrack --lat 35.58 --long 139.48 track 25544 --interval 5 --duration 300
```

#### `passes` - Predict satellite passes

Predicts upcoming satellite passes using local TLE data. This command calculates passes locally without requiring an API key.

Options:

- `--days <n>` - Number of days to search for passes (default 1)
- `--elev <n>` - Filters out passes whose maximum elevation is below this value, in degrees (default 10)
- `--horizon <n>` - Passes start and end when the satellite crosses this elevation, in degrees (default 0)

```sh
sattrack --lat 35.58 --long 139.48 passes 25544
sattrack --lat 35.58 --long 139.48 passes 25544 --days 2 --elev 20
sattrack --lat 35.58 --long 139.48 passes 25544 --horizon 10  # useful for antenna tracking
```

### Global Options

- `--config <path>` - Path to configuration file (default: ~/.sattrack.toml)
- `--lat <lat>` - Latitude of ground station in decimal degrees
- `--long <long>` - Longitude of ground station in decimal degrees
- `--alt <alt>` - Altitude above sea level in meters
- `-v, --verbose` - Display debugging information

## Dependencies

There was no need for me to reinvent the wheel, so I used existing libraries wherever it made sense to do so.

SatTrack depends on these libraries:

1. [CMake](https://cmake.org/) is used for the build system.
2. [RapidJSON](https://rapidjson.org) is used to parse JSON.
3. [curlpp](https://github.com/jpbarrette/curlpp) wraps [libcurl](https://curl.se/libcurl/) in a C++ interface for making HTTPS requests.
4. [CLI11](https://github.com/CLIUtils/CLI11) for command-line argument parsing.
5. [spdlog](https://github.com/gabime/spdlog) for fast, header-only logging.
6. [fmt](https://github.com/fmtlib/fmt) for modern C++ string formatting (required by spdlog).
7. [Asio](https://think-async.com/Asio/) for cross-platform async I/O, signals, timers, and serial ports (standalone, non-Boost version).
8. [Google Test](https://github.com/google/googletest) for unit testing (only needed for running tests).

SatTrack also depends on Howard Hinnant's [date library](https://howardhinnant.github.io/date/date.html) to support C++20 date/time functions that are not yet well supported by most compilers.

Here is the copyright statement for the date library:

```cpp
// The MIT License (MIT)
//
// Copyright (c) 2015, 2016, 2017 Howard Hinnant
// Copyright (c) 2016 Adrian Colomitchi
// Copyright (c) 2017 Florian Dang
// Copyright (c) 2017 Paul Thompson
// Copyright (c) 2018, 2019 Tomasz Kamiński
// Copyright (c) 2019 Jiangang Zhuang
```

## License

MIT License

Copyright 2025, Andrew C. Young <andrew@vaelen.org>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

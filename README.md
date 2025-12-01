# SatTrack

SatTrack is a command-line utility for tracking satellites and predicting when they will be visible from a given location. It's targeted mainly at the amateur radio and amateur satellite communities.

## Build instructions

Install dependencies:

```sh
sudo apt install build-essential cmake libcurl4-openssl-dev libcurlpp-dev rapidjson-dev libcli11-dev libgtest-dev
```

Build the application:

```sh
mkdir build
cd build
cmake ..
make
```

## Implementation details

SatTrack is a command-line application built using C++20. I wrote it largely as an opportunity to better familiarize myself with the features of C++11, C++14, C++17, and C++20 since I originally learned C++ back in the late 90s.

I used Claude Code to help speed up my implementation of SatTrack by helping me maintain my build scripts and perform research on C++20 features that I wasn't aware of, but the code was written by hand.

## Data sources

SatTrack gets both TLE and pass data for the satellites it is tracking from [N2YO](https://www.n2yo.com). You will need an API key to use the application, which you can get by signing up for a free account at N2YO.

I was originally going to implement my own pass calculator, but I haven't done that yet since this is just a small toy project at the moment. I'll try to do it later on though.

## Configuration

All configuration options can be provided on the command line, but it will be easier for you if you create a config file to store common options, like your API key and the list of satellites you care about. The config file is stored in TOML format and the options available are the same as the command line options. The config file is located at ~/.sattrack.toml by default, but you can also supply a path to the config file using the `--config` command line option.

Here is an example config file:
```
key = 'YOUR-API-KEY-HERE'
lat = 35.58
long = 139.48
id = [25544,36122,28895,27848,27844,32785,22825,7530,25397,24278]
alt = 30
```

## Usage

### Commands

SatTrack provides the following commands:

#### `passes` - Display future satellite passes

Shows upcoming passes for all configured satellites, sorted by start time.

Options:

- `--days <n>` - Number of days worth of passes to display (default 1, max 10)
- `--elev <n>` - Minimum pass elevation above the horizon in degrees (default 10, max 90)

```sh
sattrack passes

Loading Passes............. Done.

Upcoming Passes:
        Satellite                   Start                      End               Duration        Starts In        Start Az      End Az       Max Az     Max Elev
------------------------- ------------------------- ------------------------- -------------- ------------------ ------------ ------------ ------------ ----------
 27844 CUTE-1              2025-11-21 16:19:25 UTC   2025-11-21 16:33:15 UTC     13m 50s      0h  4m 28s        135.66 SE      1.09 N      67.50 ENE    27.12
 28895 CUBESAT XI-V        2025-11-21 16:42:45 UTC   2025-11-21 16:53:15 UTC     10m 30s      0h 27m 48s        348.77 N     241.77 WSW   294.94 WNW    12.84
```

#### `tle` - Display TLE data

Displays the Two-Line Element set for configured satellites.

```sh
sattrack tle
```

#### `tle parse` - Parse and display orbital elements

Parses the TLE and displays the orbital elements in a human-readable format.

```sh
sattrack tle parse
# Output:
# ISS (ZARYA)
#   NORAD ID: 25544
#   Classification: U
#   Designator: 98067A
#   Epoch: 2025-11-21 12:34:56 UTC
#   Inclination: 51.6414 deg
#   Right Ascension of Ascending Node: 247.4825 deg
#   Eccentricity: 0.0006102
#   Argument of Perigee: 52.7634 deg
#   Mean Anomaly: 64.3421 deg
#   Mean Motion: 15.50094728 revs per day
#   ...
```

#### `geo` - Get geodetic location

Calculates the geodetic location (latitude, longitude, altitude) of a satellite at a given time.

```sh
sattrack geo 25544 "2025-12-01 12:00:00"
# Output:
# Satellite: ISS (ZARYA)
#   Latitude: 23.45 deg
#   Longitude: -78.12 deg
#   Altitude: 420.5 km
```

#### `add` - Add a satellite to track

Adds a satellite by its NORAD ID to your configuration file.

```sh
sattrack add 25544
```

#### `remove` - Remove a satellite

Removes a satellite from your configuration file.

```sh
sattrack remove 25544
```

#### `key` - Set API key

Sets your N2YO API key in the configuration file.

```sh
sattrack key YOUR-API-KEY-HERE
```

### Global Options

- `--config <path>` - Path to configuration file (default: ~/.sattrack.toml)
- `--key <key>` - N2YO API key
- `--id <id>` - NORAD ID of satellite to track (can be specified multiple times)
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

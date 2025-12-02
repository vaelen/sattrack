# SatTrack

SatTrack is a command-line utility for tracking satellites and predicting when they will be visible from a given location. It's targeted mainly at the amateur radio and amateur satellite communities.

It has the following features:

1. Download TLE data from the Internet
2. Display orbital elements for one or more satellites
3. Calculate the position of a satellite at a given time
4. Displaay upcomming satellite passes.
   - The built-in pass calculator is still a work in progress, but the app can download pass information from the Internet.

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
make test
```

## Implementation details

SatTrack is a command-line application built using C++20. I wrote it largely as an opportunity to better familiarize myself with the features of C++11, C++14, C++17, and C++20 since I originally learned C++ back in the late 90s.

I used Claude Code to help speed up my implementation of SatTrack by helping me maintain my build scripts, write test cases, and perform research on C++20 features that I wasn't aware of, but the actual code was written by me.

## Data sources

SatTrack uses [CelesTrak](https://celestrak.org/) as its primary source for TLE data. CelesTrak is a non-profit service that provides satellite information for free. If you feel so inclined, they take donations to help cover their costs.

In addition, SatTrack can get both TLE and pass data from [N2YO](https://www.n2yo.com). If you want to use the N2YO APIs you will need an API key, which you can get by signing up for a free account at N2YO.

## Configuration

All configuration options can be provided on the command line, but you can also create a config file to store common options, like your N2YO API key and the list of satellites you care about. The config file is stored in TOML format and the options available are the same as the command line options. The config file is located at `~/.sattrack.toml` by default, but you can also supply a path to the config file using the `--config` command line option.

Here is an example config file:

```toml
# Station Location
lat = 35.58
long = 139.48
alt = 30

info.id = 25544
tle.id = 25544

# N2YO Settings
n2yo.key = 'YOUR-API-KEY-HERE'
n2yo.passes.id = [25544,36122,28895,27848,27844,32785,22825,7530,25397,24278]
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

# This doesn't work yet
# $ sattrack passes 25544

# Get the list of upcoming passes from N2YO
$ sattrack n2yo passes 25544

Loading Passes.... Done.

Upcoming Passes:
        Satellite                   Start                      End               Duration        Starts In        Start Az      End Az       Max Az     Max Elev 
------------------------- ------------------------- ------------------------- -------------- ------------------ ------------ ------------ ------------ ----------
 25544 SPACE STATION       2025-12-02 08:54:05 UTC   2025-12-02 09:04:40 UTC      10m 35s        0h 18m 36s      305.89 NW    146.38 SE    227.42 SW     44.89   
 25544 SPACE STATION       2025-12-02 23:58:40 UTC   2025-12-03 00:08:45 UTC      10m  5s       15h 23m 11s      200.24 SSW    61.04 ENE   129.79 SE     24.05   
 25544 SPACE STATION       2025-12-03 01:35:15 UTC   2025-12-03 01:45:35 UTC      10m 20s       16h 59m 46s      250.99 WSW    39.72 NE    324.91 NW     29.10   
 25544 SPACE STATION       2025-12-03 06:29:55 UTC   2025-12-03 06:39:10 UTC       9m 15s       21h 54m 26s      325.74 NW     84.77 E      25.10 NNE    13.44   
 25544 SPACE STATION       2025-12-03 08:06:25 UTC   2025-12-03 08:17:05 UTC      10m 40s       23h 30m 56s      311.42 NW    133.40 SE    179.33 S      88.33   
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

#### `n2yo` - N2YO API commands

Commands that fetch data from the N2YO API. These require an API key.

##### `n2yo passes` - Display future passes

Shows upcoming passes for satellites, sorted by start time.

Options:
- `--key <key>` - N2YO API key
- `--days <n>` - Number of days worth of passes to display (default 1, max 10)
- `--elev <n>` - Minimum pass elevation above the horizon in degrees (default 10, max 90)

```sh
sattrack n2yo passes --key YOUR-API-KEY 25544
```

##### `n2yo tle` - Display raw TLE data

Fetches and displays Two-Line Element sets from N2YO.

```sh
sattrack n2yo tle --key YOUR-API-KEY 25544
```

##### `n2yo info` - Display orbital elements

Fetches TLE from N2YO and displays parsed orbital elements.

```sh
sattrack n2yo info --key YOUR-API-KEY 25544
```

##### `n2yo update` - Update local TLE from N2YO

Downloads TLE data from N2YO and stores it in the local database.

```sh
sattrack n2yo update --key YOUR-API-KEY 25544 28895
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
5. [Google Test](https://github.com/google/googletest) for unit testing (only needed for running tests).

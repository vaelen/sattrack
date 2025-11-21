# SatTrack

SatTrack is a command-line utility for tracking satellites and predicting when they will be visible from a given location. It's targeted mainly at the amateur radio and amateur satellite communities.

## Build instructions

Install dependencies:

```sh
sudo apt install build-essential cmake libcurl4-openssl-dev libcurlpp-dev rapidjson-dev libcli11-dev
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

SatTrack get both TLE and pass data for the satellites it is tracking from [N2YO](https://www.n2yo.com). You will need an API key to use the application, which you can get by signing up for a free account at N2YO.

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

If you have already set up your configuration file, you can simply run:
```
$ sattrack passes

Loading Passes............. Done.

Upcoming Passes:
        Satellite                   Start                      End               Duration        Starts In        Start Az      End Az       Max Az     Max Elev 
------------------------- ------------------------- ------------------------- -------------- ------------------ ------------ ------------ ------------ ----------
 27844 CUTE-1              2025-11-21 16:19:25 JST   2025-11-21 16:33:15 JST   13 min 50 s     0 h  4 min 28 s   135.66 SE      1.09 N      67.50 ENE    27.12   
 28895 CUBESAT XI-V        2025-11-21 16:42:45 JST   2025-11-21 16:53:15 JST   10 min 30 s     0 h 27 min 48 s   348.77 N     241.77 WSW   294.94 WNW    12.84   
 27848 CUBESAT XI-IV       2025-11-21 16:56:15 JST   2025-11-21 17:11:10 JST   14 min 55 s     0 h 41 min 18 s   156.40 SSE   351.58 N      72.83 ENE    62.60   
 7530  OSCAR 7             2025-11-21 17:03:50 JST   2025-11-21 17:26:00 JST   22 min 10 s     0 h 48 min 53 s   162.41 SSE   343.10 NNW   250.54 WSW    85.23   
 27844 CUTE-1              2025-11-21 17:58:55 JST   2025-11-21 18:13:15 JST   14 min 20 s     1 h 43 min 58 s   189.89 S     335.34 NNW   261.64 W      32.14   
 27848 CUBESAT XI-IV       2025-11-21 18:37:55 JST   2025-11-21 18:50:05 JST   12 min 10 s     2 h 22 min 58 s   213.28 SW    322.43 NW    266.46 W      13.40   
 7530  OSCAR 7             2025-11-21 19:00:00 JST   2025-11-21 19:17:10 JST   17 min 10 s     2 h 45 min  3 s   215.71 SW    320.93 NW    267.04 W      14.00   
 24278 JAS 2 (FO-29)       2025-11-21 19:20:50 JST   2025-11-21 19:36:35 JST   15 min 45 s     3 h  5 min 53 s    32.93 NE    140.68 SE     86.00 E      14.52 
```

## Dependencies

There was no need for me to reinvent the wheel, so I used existing libraries wherever it made sense to do so.

SatTrack depends on these libraries:

1. [CMake](https://cmake.org/) is used for the build system.
2. [RapidJSON](https://rapidjson.org) is used to parse JSON.
3. [curlpp](https://github.com/jpbarrette/curlpp) wraps [libcurl](https://curl.se/libcurl/) in a C++ interface for making HTTPS requests.
4. [CLI11](https://github.com/CLIUtils/CLI11) for command-line argument parsing.

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
sattrack passes
```

## Dependencies

There was no need for me to reinvent the wheel, so I used existing libraries wherever it made sense to do so.

SatTrack depends on these libraries:

1. [CMake](https://cmake.org/) is used for the build system.
2. [RapidJSON](https://rapidjson.org) is used to parse JSON.
3. [curlpp](https://github.com/jpbarrette/curlpp) wraps [libcurl](https://curl.se/libcurl/) in a C++ interface for making HTTPS requests.
4. [CLI11](https://github.com/CLIUtils/CLI11) for command-line argument parsing.

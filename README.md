# SatTrack

SatTrack is a small utility for keeping track of satellites that are currently in sight of a given location. There are other more feature filled satellite trackers out there. This utility is meant to sit in your tray and notify you when satellites that you are interested in are going to be in view. It's targeted mainly at the amateur radio and amateur satellite communities.

## Build instructions

Install dependencies:

```sh
sudo apt install build-essential cmake libgtkmm-3.0-dev libcurl4-openssl-dev libcurlpp-dev rapidjson-dev libcli11-dev
```

Build the application:

```sh
mkdir build
cd build
cmake ..
make
```

## Implementation details

SatTrack is a GTK++ application and is built using C++20. I wrote it largely as an opportunity to better familiarize myself with the features of C++11, C++14, C++17, and C++20 since I originally learned C++ back in the late 90s.

I used Claude Code to help speed up my implementation of SatTrack by helping me maintain my build scripts, write tests, perform repetative tasks, and perform research, but the application itself was mostly written by hand.

## Data sources

SatTrack uses the standard TLE format for orbital elements of the satellites it is tracking which it downloads from [N2YO](https://www.n2yo.com). NY2YO provides a [RESTful API](https://www.n2yo.com/api/) for retrieving TLE data.

## Dependencies

There was no need for me to reinvent the wheel, so I used existing libraries wherever it made sense to do so.

SatTrack depends on these libraries:

1. [CMake](https://cmake.org/) is used for the build system.
2. [GTKmm](https://gtkmm.gnome.org/en/index.html) for the UI elements.
3. [RapidJSON](https://rapidjson.org) is used to parse JSON.
4. [curlpp](https://github.com/jpbarrette/curlpp) wraps [libcurl](https://curl.se/libcurl/) in a C++ interface for making HTTPS requests.
5. [CLI11](https://github.com/CLIUtils/CLI11) for command-line argument parsing.

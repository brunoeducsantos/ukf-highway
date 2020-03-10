# Unscented Kalman Filter in Highway simulation


<img src="media/ukf_highway_tracked.gif" width="700" height="400" />

## Overview
This project is an implementation of Unscented Kalman Filter to estimate the state of multiple cars on a highway using noisy lidar and radar measurements. 

The viewer scene is centered around the ego car and the coordinate system is relative to the ego car as well. The ego car is green while the  other traffic cars are blue. The traffic cars will be accelerating and altering their steering to change lanes. Each of the traffic car's has it's own UKF object generated for it, and will update each indidual one during every time step. 

The red spheres above cars represent the (x,y) lidar detection and the purple lines show the radar measurements with the velocity magnitude along the detected angle. The Z axis is not taken into account for tracking, so you are only tracking along the X/Y axis.


## Dependencies
* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
 * PCL 1.2

## Usage

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./ukf_highway`

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html) as much as possible.

## Generating Additional Data

If you'd like to generate your own radar and lidar modify the code in `highway.h` to alter the cars. Also check out `tools.cpp` to change how measurements are taken, for instance lidar markers could be the (x,y) center of bounding boxes by scanning the PCD environment and performing clustering. 


## Disclamer
This project was cloned from [Udacity UKF project](https://github.com/udacity/SFND_Unscented_Kalman_Filter) in the context of [Sensor Fusion Engineer nanodegree](https://www.udacity.com/course/sensor-fusion-engineer-nanodegree--nd313).


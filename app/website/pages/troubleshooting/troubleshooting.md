---
title: Troubleshooting
nav_order: 7
has_childre: False
---

# Troubleshooting

This page provides simple troubleshooting steps for analysis with GelBox 3. Clicking on any of the images on this page will open a larger version in a new browswer window.


Here's an example of a situation that you might encounter with your own gels. Please note that the fitted function does not properly represent the density profile. The application only captures one peak and as a result there is only one non-zero Gaussian function.

The problem here is that the software hasn't managed to find a good fit. There are two distinct possibilities.

+ There is not  a good mathematical fit to the experimentally-measured profile, and
+ There is a good fit and GelBox 3 could not find it

<a href="media/single_curve.png" target="_blank">![Single curve](media/single_curve.png)</a>

It is not always easy to identify the situation. One possible first step is to visualize the fitting process to observe how the function develops through iterations. Draw Fitting checkbox, shown in red rectangle, lets users to visualize fitting. As soon as the box is checked, fitting starts over. Please note that this might take a little bit longer as each fitting iteration is simultenously plotted.

draw_fitting.mp4

The first couple iterations show that the peak locations are not accurately predicted by the application. Click the Fitting Parameters in the Fitting panel for the starting parameter estimates. Although the peak locations are roughly around 20 and 60 for band 1 and band 2, the estimated parameters are 100 and 80, respectively.

<a href="media/change_parameters.png" target="_blank">![Change parameters](media/change_parameters.png)</a>

Change the peak locations for both band 1 and band 2 and click the Update Fitting Parameters button.

<a href="media/parameters_changed.png" target="_blank">![Parameters changed](media/parameters_changed.png)</a>

Please note that the upon changing the peak locations, GelBox 3 succesfully detected the two peaks. You can look for better fits by changing the size and position the region of interest (ROI) box.

<a href="media/drag_box_down.png" target="_blank">![Drag box down](media/drag_box_down.png)</a>

Further improvement can be made by changing the shape parameter of the second band using the Fitting Parameters.

<a href="media/final_analysis.png" target="_blank">![final analysis](media/final_analysis.png)</a>

Through these troubleshooting steps, an acceptable R-squared value is obtained for the fitted function. In addition, both bands are represented with non-zero Gaussian functions.

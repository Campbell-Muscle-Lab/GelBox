---
title: Troubleshooting
nav_order: 7
has_childre: False
---

# Troubleshooting

This page provides simple troubleshooting steps for analysis with GelBox 3. Clicking on any of the images on this page will open a larger version in a new browser window.

Here's an example of a situation you might encounter with your gels. Please note that the fitted function does not properly represent the density profile. The application only captures one peak, and as a result, there is only one non-zero Gaussian function.

The problem is that the software has not found a good fit. There are two distinct possibilities.

+ There is not  a good mathematical fit to the experimentally-measured profile, and
+ There is a good fit, and GelBox 3 could not find it

<a href="media/single_curve.png" target="_blank">![Single curve](media/single_curve.png)</a>

It is not always easy to identify the situation. One possible first step is visualizing the fitting process to observe how the function develops through iterations. Draw Fitting checkbox, shown in the red rectangle, lets users visualize fitting. As soon as the box is checked, the fitting starts over. Please note that this might take a little bit longer as each fitting iteration is simultaneously plotted.

<video src="https://github.com/Campbell-Muscle-Lab/GelBox/assets/98066302/95cf4387-1897-48d9-ac23-2b83ad690f51" controls="controls" style="max-width: 730px;"></video>

The first couple of iterations shows that the application does not accurately predict the peak locations. Click the Fitting Parameters in the Fitting panel for the starting parameter estimates. Although the peak locations are roughly around 20 and 60 for band 1 and band 2, the estimated parameters are 100 and 80, respectively.

<a href="media/change_parameters.png" target="_blank">![Change parameters](media/change_parameters.png)</a>

Change the peak locations for both band 1 and band 2 and click the Update Fitting Parameters button.

<a href="media/parameters_changed.png" target="_blank">![Parameters changed](media/parameters_changed.png)</a>

Please note that after changing the peak locations, GelBox 3 successfully detected the two peaks. You can look for better fits by changing the size and position the region of interest (ROI) box.

<a href="media/drag_box_down.png" target="_blank">![Drag box down](media/drag_box_down.png)</a>

Change the peak estimate of the peak one and constrain it to capture the amplitude accurately. Use the constrain checkbox, shown in the red rectangle.

<a href="media/constrain_parameter.png" target="_blank">![Drag box down](media/constrain_parameter.png)</a>

After constraining the peak amplitude, the new fit accurately predicts the peak amplitude.

<a href="media/constrained_parameter_new_fit.png" target="_blank">![Drag box down](media/constrained_parameter_new_fit.png)</a>

Further improvement can be made by changing the shape parameter of the second band using the Fitting Parameters.

<a href="media/final_analysis.png" target="_blank">![final analysis](media/final_analysis.png)</a>

These troubleshooting steps obtain an acceptable R-squared value for the fitted function. In addition, both bands are represented with non-zero Gaussian functions.
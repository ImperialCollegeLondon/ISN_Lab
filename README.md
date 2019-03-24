# InterventionalSystemsNeuroscience
The aim for this branch is to develop and write MATLAB functions that perform common processes in most of our analysis scripts.
Anything that is often repeated can be turned into a function. This is not only beneficial in the sense that several lines of code can become one line in our scripts, it also means that, especially fr simple functions, we can obsfuscate the code allowing MATLAB to run faster.

# Development Process
Write MATLAB functions '*function*.m' that are as robust and universal as possible.

Make sure the descriptions are accurate and clear enough to allow anybody to use the function.

Once the '.m' function has been written we can save it in the 'DevFunctions' folder, this allows edits to be made if necessary.

If you are happy with the function, you can obsfuscate it in MATLAB by running pcode('*function*').

The resulting obsfucated function '*function*.p' can be uploaded into the 'Functions' folder.

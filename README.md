Hello there!

I provide here a quick sumary for you to be able to fully understand the architecture of the SaaS if you're not friendly with web developpement and even, add some features to the existing project.

The web application is handle by the flask micro-framework, available as a python libraries which provides, the direction for the index.html and handle the HTTPs requests send by the fetch API method define in Js.

Some particularities of the formats conversion is in the mol file treatment. Indeed, the mol file is defined as an text which is composed of several pre-informations (1st : chemical formula, 
2d: infos on the way of obtention of the file, 3th: an empty line which compose the "header lines", 4th: counts line, 5th and more: beginning of atom block and then, the bond block after)

Exemple that works: 

C3H7O
     RDKit          3D

 12 11  0  0  0  0  0  0  0  0999 V2000
    1.0557    0.5879    1.4293 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8945   -0.1960    0.1599 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2889    0.7287   -0.8988 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1474    1.8082   -1.1006 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7799    0.0074    2.3214 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1091    0.9087    1.5244 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4037    1.5098    1.3926 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8773   -0.5387   -0.2078 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.2380   -1.0641    0.2429 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6635    1.1165   -0.4461 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0675    0.1766   -1.8384 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7485    2.6428   -0.7437 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  2  8  1  0
  2  9  1  0
  3 10  1  0
  3 11  1  0
  4 12  1  0
  §§§§
  M  END
  §§§§

  Exemple that will bring an error:

     RDKit          3D

 12 11  0  0  0  0  0  0  0  0999 V2000
    1.0557    0.5879    1.4293 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8945   -0.1960    0.1599 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2889    0.7287   -0.8988 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1474    1.8082   -1.1006 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7799    0.0074    2.3214 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1091    0.9087    1.5244 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4037    1.5098    1.3926 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8773   -0.5387   -0.2078 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.2380   -1.0641    0.2429 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6635    1.1165   -0.4461 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0675    0.1766   -1.8384 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7485    2.6428   -0.7437 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  2  8  1  0
  2  9  1  0
  3 10  1  0
  3 11  1  0
  4 12  1  0
  §§§§
  M  END
  §§§§

Some informations about why we need to keep this form.
The python code have 4 steps of making sure the mol file is fonctionable and the verification with the molecular formula in regards of the graph is one (cf VERIFICATION in class_molecule).

##########################
####Modification parts####
##########################

To add new features, the simplest way is to first code in python in a separate file what you want to do in back-end after the input of the user.
Then when that is ready, you need to recuperate the input users in the dict_inputs define in the @app_route.
For the js script, if you did what is before, the results should appear in the JSON dict of outputs which is contain in "data".
You then just need to assign a js constant to this value and make it print on the front-end.
Finally, in the HTML script, you need to include a input of desire type (text, checkbox, boutton, ...) in the form for the users to be able to inputs some values.



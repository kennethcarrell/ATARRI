---
title: 'ATARRI: A TESS Archive RR Lyrae Classifier'
tags:
- Python
- astronomy
- TESS
- GUI
- RR Lyrae
authors:
- name: Kenneth W. Carrell
  orcid: 0000-0002-6307-992X
  affiliation: 1
affiliations:
- name: Angelo State University
  index: 1
date: 06 May 2021
bibliography: paper.bib

aas-doi: xxxx/xxxx
aas-journal: The Astrophysical Journal Letters
---

# Summary

Radially pulsating variable stars change in brightness over time. For
RR Lyrae type stars specifically, this change occurs over
approximately half a day. Surveys collect data on these stars
over long periods of time and can be used to quantify things such as
the exact period and amplitude of the change in brightness. Automated
procedures are normally used to analyze large datasets of known variable
stars, and to look for other potentially variable objects. However,
changes in the period and/or amplitude of the pulsation cycle require
a more detailed analysis.

The ATARRI GUI is designed as a user-friendly interface to display
information about known or suspected RR Lyrae variable stars. The data
comes from the Transiting Exoplanet Survey Satellite [TESS, @tess] data
archive, which upon its completion will have 
imaged almost the entire sky twice. Being an all-sky survey, it is 
obviously very large in terms of types of objects observed as well as
amount of data stored. Downloading data is accomplished with the
`search_tesscut` functionality [@tesscut] available in the
`lightkurve` python package [@lightkurve], which also has methods used 
for the analysis presented in the GUI.

Non-professional astronomers, and specifically
undergraduate students in STEM fields, can be easily trained to
inspect the information in the GUI and provide some basic information
about each object.

# Statement of need

Astronomical surveys produce very large datasets and tools must be
developed to process and analyze this information. Many of these
tools, if not most, assume that the behavior of individual objects are
so similar that a single process or routine can describe them. While
this may be true, and characterizing global properties of certain
objects is important to understanding them, things that exhibit
strange behavior are also very interesting and often overlooked in
large surveys. Furthermore, as useful as algorithms are at work of
this nature, including machine learning techniques, training an
algorithm to identify odd behavior or subtle changes ranges from very
difficult to impossible. Ideally, a professional astronomer would
inspect every object in the dataset to verify data quality and
analysis, and thereby provide a check on computed values and
classification. In practice, however, datasets are much too large for
even a team of astronomers to individually inspect every single
object.

One solution to this is to train non-professionals to perform specific
tasks so that a human looks at the data and can check certain,
specific things. So-called citizen science projects have been an
enormous help in some areas of research to quickly and efficiently
analyze large datasets. With the ATARRI software, we have taken a
similar approach to classify and analyze a specific type of variable
star.

# Acknowledgements

Funding was provided to KC through a grant from the Angelo State
University Faculty Research Enhancement Program.

# References

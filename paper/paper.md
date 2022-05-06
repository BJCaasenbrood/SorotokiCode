---
title: 'SOROTOKI: A Soft Robotics Toolkit for Matlab '
tags:
  - Matlab
  - soft robotics
  - finite elements
  - continious dynamics
  - topology optimization
authors:
  - name: Brandon J. Caasenbrood
    orcid: 0000-0000-0000-0000
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Alexander Y. Pogromsky
    orcid: 0000-0000-0000-0000
    affiliation: "1,2"
  - name: Henk Nijmeijer
    affiliation: "1"
affiliations:
 - name: Lyman Spitzer, Jr. Fellow, Princeton University
   index: 1
 - name: Institution Name
   index: 2
date: 23 April 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
In the past years, the field of soft robotics has spread rapidly among the scientific community. Unlike rigid robotics, soft robots are purposefully composed of soft materials allowing for intrinsic compliance and safety, a rich family of continuum-bodied motion, and environmental resilience. The main inspiration for soft robots stems from biology with the aim to achieve similar performance and dexterity as biological creatures. Although the field has made major steps have been towards bridging biology and robotics, its innate infinite-dimensionality poses significant challenges on modeling and control. In terms of performance, soft robots are easily outclassed by their rigid counterparts and consequently lack the transferability to industry [5]. The diligence of achieving similar precision and speed to nowadays’ rigid robots, and ultimately nature, stresses the paramount importance of novel modeling and control strategies tailored for soft robotics. 

`SOROTOKI` is an open-source MATLAB toolkit for soft robotics that includes an array of tools for design, modeling, and control. Due to its scientific diversity, it can be challenging for researchers to quickly familiarize themselves with multiple scientific disciplines. With the aim to lower the threshold, Sorotoki aims to incorporate multiple layers of soft robotics research into one toolkit. Examples include: continuum mechanics, dynamic systems and control theory, topology optimization, computer graphics, and much more to come! The combination provides a highly flexible modeling environment and hopefully aids the development of novel soft robotic research.

# Mathematics
Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures
Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }
# Acknowledgements

# References
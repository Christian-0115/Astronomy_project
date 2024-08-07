
# Project Title

This is an astronomy project that I started on May 2022 under the mentorship of Dr. Gregory Simonian at Bradley University. The research explores the nature of binary star systems, blending astronomy with the power of computer science. My research was focused on the orbital properties of binary systems.

The Kepler satellite found a population of 27 stars rotating more quickly than expected for single stars.  After an initial analysis of the rapidly rotating stars with more detailed observations was inconclusive, we decided to build a more sophisticated model in Python that we can compare the data to support or reject our binary hypothesis. Our hypothesis is that these are spectroscopic binary systems, where tidal interactions between the stars keep them rotating quickly. Synchronous tidal locking occurs when a star's orbital period is the same as its rotational period.

Kepler’s third law of planetary motion allows astronomers to derive the measurement of the mass of both stars. Using Kepler’s 3rd law of planetary motion, our Python model produced radial velocity curves that model the behavior of a binary star system based on known quantities about the stars as well as several assumptions. 

The known quantities of our target stars include the rotational period acquired from McQuillan et al (2014) and the mass acquired from Berger et al (2020). Our Python model used a statistics Chi-square test to  find the best-fit radial velocity curve that best fits the movement of a binary star system. We then compare our modeled radial velocity curve with the actual measured radial velocity times, measured by the Modspec instrument at the Michigan-Dartmouth-MIT (MDM) observatory.


## Authors

- [@gregosimo](https://www.github.com/gregosimo)
- [@Christian-0115](https://wwwgithub.com/Christian-0115)

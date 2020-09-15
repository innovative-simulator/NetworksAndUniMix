# NetworksAndUniMix
Compare disease transmission in a social network with universal mixing

NetworksAndUniMix.nlogo is a computer simulation written in NetLogo. It simulates the transmission of an infectious disease. In common with most mathematical epidemiologists' models (e.g. the SIR model of Kermack & McKendrick 1927), it divides a population into compartments according to their disease state: Susceptible, Exposed (but not yet infectious), Infectious, or Removed from being able to transmit it (because you are Recovered and immune, or Dead, or just isolating). Susceptible people are infected through social contacts with infectious people. The rate at which this occurs depends on their social contact rate (how many social contacts per time step they have), and on the proportion of those contacts that are with infectious people. A common assumption in the mathematical models is that of universal mixing: that your social contacts are with anyone in the population, irrespective of where they are in space or whether you have met them before. A common criticism of these models, therefore is that real people live in spatial and social networks, social contacts are limited to some kind of network's neighbours, and disease transmission occurs via contacts with these neighbours, not with just anyone. This NetLogo program compares transmission via universal mixing with transmission via a social network, while keeping other parameters (especially the contact rate) constant.

There are a variety of social network structures to choose from - common library models, such as Strogatz & Watts small world, Barabasi-Albert scale free, Erdos-Renyi random, and regular networks such as 2-D grids and rings. We also include Hamill & Gilbert's (2009) social circles network, and a network based on real contact data published by Kissler et al. (and collected in conjunction with the BBC's 2018 TV programme, Pandemic). 

For the Kissler network, you will need to download the two .csv datafiles and place them in the same folder as NetworksAndUniMix.nlogo.

The .xlsx file contains the results of experiments run using BehaviorSpace.

The .pdf file are the slides from my presentation to Social Simulation Week 2020, given online on 16th September 2020.

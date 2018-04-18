# Introduced_Tradeoffs
Code for reproducing figures in "Introduced species that overcome life history tradeoffs can cause native extinctions"

This repository contains three functions.

1. GeneratePersistentCommunities.m
This function creates metacommunities with a given number of stably coexisting species. It randomly generates parameters for a Tilman-1994-style metacommunity of multiple species interacting according to a strict competitive hierarchy for space, and then checks to see if they can stable coexist. It retains those communities where all species persist.

2. ForwardSimulate
This function models the time-evolution of the communities described above.

3. Fig_Timeseries_of_extinctions
This function produces Figure 1C from the manuscript. It uses the persisting communities identified by GeneratePersistentCommunities, and solves for their abundances through time using ForwardSimulate.

Adding a new field, mp0, to particle data, which represents the initial mass
of a star particle with a flag -DINIT_STELLAR_MASS
When you activate -DSTELLAR_POPULATION_MASS, a star particle with mass smaller than 100 Msun will assign initial stellar population mass, msp0. Luminosity and stellar winds energy estimated based on this mass, while mass-loss is calculated based on mp0.


# Hanging films

Numerical method implementation to study films hanging on the underside of an inclined plate.
Both a Benney and weighted residual integral boundary layer model are used and compared.
This work is done in collaboration with Ruben Tomlin and Demetrios Papageorgiou.

### Acceptance tests.

The coding solution should:

- apply a numerical scheme to the given model problem
- save the resulting output
- be flexible in which model is being used
- be flexible in which numerical scheme is being used

### Unit tests.

Running the unit tests in Matlab,

    runtests;

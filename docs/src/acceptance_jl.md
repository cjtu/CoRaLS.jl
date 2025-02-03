# acceptance.jl: The Main Monte Carlo Loop

The main `acceptance()` loop computes many cosmic ray impacts to a supplied region of the lunar surface and observed with a chosen spacecraft and triggering method, returning an `Acceptance` object containing the number of detected events and final acceptance in [$km^2$ sr] units.

```@autodocs
Modules = [CoRaLS]
Pages = ["acceptance.jl"]
```

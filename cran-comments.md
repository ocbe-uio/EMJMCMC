Dear CRAN Team,

this is a resubmission of package 'EMJMCMC' (not yet on CRAN). We have addressed all points raised by the CRAN team in the previous submission, except for the replacement of the <<- operators, which we understand are an essential part of the ReferenceClasses system. From ?ReferenceClasses:

"Fields may be modified in a method by using the non-local assignment operator, <<-, as in the $edit and $undo methods in the example below. Note that non-local assignment is required: a local assignment with the <- operator just creates a local object in the function call, as it would in any R function. When methods are installed, a heuristic check is made for local assignments to field names and a warning issued if any are detected."

We have made attempts to replace this operator, but it is not possible to do so without breaking the functionality of the package. We understand that <<- is the only way to modify class fields in ReferenceClasses, but look forward to hearing from you if there is a way to address this issue that doesn't involve completely replacing class systems (which sounds unnecessary, given ReferenceClasses is a standard part of R).

Best, Waldir

# Package EMJMCMC 1.5.0

Reporting is done by packager version 1.15.2


## Test environments
- R version 4.3.1 (2023-06-16)
   Platform: x86_64-pc-linux-gnu (64-bit)
   Running under: Ubuntu 23.10
   ERROR: No check log found!
- win-builder (devel)

## Local test results

## Local meta results

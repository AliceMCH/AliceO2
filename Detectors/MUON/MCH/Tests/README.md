ReadMe - MCH/Tests/

There are 3 tests that run the preclustering/clustering procedure.



The first one is ClusterizerTest.
It takes as input a binary file of digits, rune the preclustering and the clustering, following the method chosen in the code ClusterizerText.cxx
There are 4 methode of clustering currently implemented:
    - Centor of Gravity
    - Mathieson Fit (called SimpleFit)
    - Gaussian Fit
    - Double Gaussian Fit
    
Usually, we run ClusterizerTest on the binary file digits.in from Laurent's CERNBox,  https://cernbox.cern.ch/index.php/s/6J4kTkmal0X2sBK

    o2-mch-clusterizer-test digits.in

This should just prompt the preclustering/clustering process, without anything fancier like graphs or anything. This test, mainly used up until end of 2019 was made to make sure that teh clustering algorithms were running without failing.
It can also be used to check that the preclusters that are prompted on console match the theoretical ones from Laurent in preclusters.out (on CERNBox as well).



The second test historically was ClusterizerValidation.
To deal with the fact that we did not have "theoretical clusters" associated with our binary digit files, we had to write some code that would generate ideal clusters, get their digits and use them as input for the preclustering/clustering (that part being similar to ClusterizerTest).

The code generates ideal events on DE819. Those events are hits in the detector leaving a Mathieson distribution of charge. You can choose to generate 1 hit or 2 hits for a single event (2 hits will be useful when trying to reconstruct non-singular clusters)

The ClusterizerValidation.cxx algorithm will generate a given number of events following the constraints on the position or charge of clusters you enter. Those events will then be reconstructed using The function TestClustering() in Validation.cxx which is essentially ClusterizerTest.cxx so you can choose there the clusterization algorithm to use.

Running ClusterizerValidation does not need an input file, just to:
    -   tweak the values in ClusterizerValidation.cxx following the number of events you want and the use you make of them
    -   choose the clustering algorithm you need in Validation.cxx (function TestClustering())
    
    o2-mch-clusterizer-validation
    
There are multiple ways you can use ClusterizerValidation:
    -   As a way to check if your clustering method leaves a bias on the residuals with respect to their position on a pad
    -   As a way to look at the distribution of the residuals over a given number of events
    -   As a way to look at the distribution of residuals (particularly the width) with respect to a fixed cluster charge
    
The input needed for these 3 use cases is detailled in ClusterizerValidation.cxx



The third and last test is ClusterizerTestTB.

ClusterizerTestTB is a direct adaptation of ClusterizerTest. The only difference is that it is used to look at binary digits files obtained from Test Beam data.
In this case, it means that contrary to ClusterizerTest using the digits.in from Laurent's CERNBox, we have an idea of where we should find our clusters thanks to the positions taken by the SiTracker.
So the difference with ClusterizerTest is that we can directly plot the residuals distribution for the TB data and see also the residuals with respect to the hit coordinate y and see if there is any bias.

The file used is a file made by Andrea called [digits-TB-run1036.bin](https://cernbox.cern.ch/index.php/s/N5jOVO9sphIQRnW)

    o2-mch-clusterizer-test-tb digits-TB-run1036.bin

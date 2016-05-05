#Analysis procedure

1. pp process reweighting

    ```
	 cd processweighting
	 rb flavorProcess.C+
	 rl checknorms.C
    ```
    use the values of beta and gamma in

    ```
	 rb skiplightreweight.C+
	 rb skiplightreweightdraw.C+
    ```
	Copy those values into `helpers/physics.h` : `vector<float> processWeights`.

2. Estimate hydjet amount.

	```
	cd mistag
	rb hydjetestimation.C+
	```
	Copy output values to `helpers/physics.h` : `vector<float> bkgfractionInNearSide `.
	Check results with
	
	```
	rb hydjetclosure.C+
	```
	
3. Get the result plots

	```
	cd asymmetry
	rb tellmetruth.C+
	```
	
	Or to test all possibilities
	
	```
	sh runall.sh
	```
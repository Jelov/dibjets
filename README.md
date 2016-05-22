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
	
	
##Remarks
* To prevent confusion, centrality bin number in the code is always 0-200 (except tagging efficiency correction functions), and shown as 0-100% only in the plotting.
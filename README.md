# Statistically Classifying Tandem Repeats Within the Human Genome
Project lead: Gary Benson\
Primary contributor: Michael Korovkin

***
#### Project Prompts
1. How different must "critical ratios" (read count inside array divided by read
count on flanks) be in order to statistically conclude whether the tandem repeat's
overall array length changed across two samples?
2. Can the length of a given tandem repeat array be predicted from ratio data?
3. Can the copy number of a tandem repeat be determined from a ratio and a pattern
length?
4. Can the genotype of a tandem repeat, given its respective human genome data, be
identified?

***
#### File Summary
```graph_color2.py```\
File useful for exploratory analysis and data visualization. Uses Plotly to display
data and allow easy manipulation and scrolling. Displays many different aspects of
the data.

```graph_color2_driver```\
Driver for <i>graph_color2.py</i>.

```simulator6.py```\
Simulates a read-mapping experiment on the human genome, using statistics emulative
of actual human genome data. Provides graphical representations of the results.

***
#### Exploratory & Preliminary Work
1. A simulator was built to emulate the experimental sequencing of the human genome
in order to obtain a simulated result of "critical mapping ratios" (calculated by
dividing the read count inside a given array region by the read count on its two
respective flanks)
2. This simulator was used to obtain expected data for the 0/0 genotype of an
arbitrary tandem-repeat-linked gene, given the expected mapping data of the gene
3. the expected data was used to calculate expected standard deviations of the
critical ratios observed at each gene array length, as well as the expected value
of the critical ratio
4. Ideal linear regression models for 0/0 genotypes were calculated from the
expected data
5. A statistical method of prediction was built to predict the read length
associated with a certain human genome dataset; the algorithm functioned with
approximately 44%

***
#### Current Work
<i>Notable conditions of data: 1) flanks are kept constant, equal to their mean
across the 0/0 genotype dataset of a simulated dataset, 2) </i>

***
#### Results & Points
1. More data would be helpful, since (as seen on the simulator) simple standard
deviation-based classification of points is not robust; there is too much noise
in the data to classify things using standard statistical methods

***
#### Conclusions
1. Problem 1 was solved through the aforementioned simulation in
```simulator6.py```
    * In summary, it was found that ratios of the same read length, corresponding
    to the same tandem repeat, must differ by more than <i>twice</i> the quantity
    ```(1.96 * standard deviation)```, where the standard deviation is obtained
    from the simulations executed in ```simulator6.py```
2. Problem 2 is relatively simple to solve if flank intersection counts for a given
tandem repeat are known; in this case, an elementary algebraic equation must be
solved in order to predict the array length of the repeat
    * If flanks are not known, but the genome/area coverage is known, then a linear
    regression model can be used to predict the array length since the flank counts
    can be approximated using the genome coverage
3. Problem 3 was straightforward to solve, since the pattern length is directly
related to the copy number of the tandem repeat
    * Knowing the pattern length and the ratio of the tandem repeat can yield an
    accurate estimation of the tandem repeat's copy number
    * if the pattern length is not known, then estimating the copy number requires
    knowledge of the given tandem repeat's genotype
    
***
#### Future Developments
1. Problem 4 will be addressed using a set of classifiers (SVM, Bayesian decision
tree, random forest); the performance and consistency of each one will be
investigated
    * The data will also be screened under statistical tests such as Chi-squared
    and F, in hopes of capturing differences in overall dataset variance
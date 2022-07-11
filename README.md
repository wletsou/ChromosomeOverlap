<h1>Instructions</h1>

This document outlines how to perform Chromosome Overlap on a single machine in a bash shell.  Users may parallelize computations on an HPC environment by selecting a range of steps to perform in each <kbd>_sub</kbd> operation.  The goal is to obtain a list of closed haplotype patterns from (filtered meta-)chromosomes from a population.  Users may download the scripts or clone them with git to a convenient location called <kbd>HOME_DIR</kbd> type the full path when needed or add them to their <kbd>PATH</kbd> variable.

<h2>Initiation</h2>

<p>Your input should be a haplotype file (<kbd>.txt</kbd>) of the following form:</p>

<table>
  <tr>
    <th>sid</th>
    <th>rs0000001_A</th>
    <th>rs0000002_G</th>
    <th>rs0000003_C</th>
    <th>...</th>
  </tr>
  <tr>
    <th>0001</th>
    <td>0</td>
    <td>0</td>
    <td>0</td>
    <td>...</td>
  </tr>
  <tr>
    <th>0001</th>
    <td>0</td>
    <td>1</td>
    <td>0</td>
    <td>...</td>
  </tr>
  <tr>
    <th>0002</th>
    <td>1</td>
    <td>0</td>
    <td>0</td>
    <td>...</td>
  </tr>
  <tr>
    <th>0002</th>
    <td>1</td>
    <td>0</td>
    <td>0</td>
    <td>...</td>
  </tr>
 </table>
 
 <p>Each row is chromosome from subject sid; note that each has two chromosomes with the same id.  Alternate alleles of SNP </kbd>rsid_Alt</kbd> are indicated with a 1 and reference alleles with a 0.  Transpose the file with one of the following commands.  Note that the sid row is purposely excluded.</p>
 
 <pre>
 <code>
 n_fields=$(awk 'END{print NF}' file)
 for i in `seq 2 $n_fields`
 do
    cut -f$i file | paste -s # paste column i as a row
 done > transpose_file
 </code>
 </pre>
 
 The output will be something like:
 
 <table>
  <tr>
    <th>rs0000001_A</th>
    <td>0</td>
    <td>0</td>
    <td>1</td>
    <td>1</td>
    <td>...</td>
  </tr>
  <tr>
    <th>rs0000002_G</th>
    <td>0</td>
    <td>1</td>
    <td>0</td>
    <td>0</td>
    <td>...</td>
  </tr>
  <tr>
    <th>rs0000003_C</th>
    <td>0</td>
    <td>0</td>
    <td>0</td>
    <td>0</td>
    <td>...</td>
  </tr>
 </table>
 
 <p>The initiation step overlaps columns of <kbd>transpose_file</kbd>.  For users not running on an HPC, the <kbd>ChromosomeOverlap_initiation_sub.sh</kbd> program loops through all <i>variable</i> chromosomes and overlaps them with each column <var>j</var> of the transposed matrix such that <var>j>i<sub>&sigma;-1</sub>.</p>
  
  <pre>
  <code>
  sh HOME_DIR/ChromosomeOverlap_initiation_sub.sh transpose_file 1 "NAME" "" "DIRECTORY" "HOME_DIR"
  </code>
  </pre>
  
  <p>The first input is the tranposed matrix, the second is <var>&sigma;-1</var> (i.e., the number of fixed chromosomes in each overlap) and <kbd>NAME</kbd> is an optional input for naming the output file; a good name is the chromosome regions, e.g., <kbd>chr11.69231642-69431642</kbd>.  The fourth empty input is for indicating a range of the <var>2N</var>-choose-<var>(&sigma;-1)</var> combinatins of fixed chromosomes. You may want to suppress output with the <kbd>> /dev/null</kbd> redirect.  In addition, if <kbd>HOME_DIR</kbd> is on your <kbd>PATH</kbd>, you do not need to specify it on the command line.</p>
  
  <p>If a range of tuples, say 0 to 49, is to be computed, the script</p> 
  
  <pre>
  <code>
  sh HOME_DIR/index2combo2.sh 0 100 1
  </code>
  </pre> 
  
  finds the tuple of <var>&sigma;-1=1</var> fixed chromosomes corresponding to the linear index 0 (i.e., the 0<sup>th</sup> tuple) when there are <var>n=100</var> total chromosomes.  This is a trivial computation when <var>&sigma;-1=1</var>, but the program easily extends to larger values of <var>&sigma;</var>.  Use of this computation is critical for parallelization.  (Note that the R version of <kbd>index2combo</kbd> is much faster; the sh version is included for portability.)</p>
  
  <p>The intermediate files generated by <kbd>ChromosomeOverlap_initiation.sh</kbd> look like </p>
  
  <table>
    <tr>
      <td>10</td>
      <td>00</td>
      <td>01,02,04,07</td>
    </tr>
    <tr>
      <td>10</td>
      <td>10</td>
      <td>03,05</td>
    </tr>
    <tr>
      <td>10</td>
      <td>01</td>
      <td>06,08</td>
    </tr>
    <tr>
      <td>10</td>
      <td>11</td>
      <td>09</td>
    </tr>
   </table>
    
<p>The first column is the number of the variable chromosome (i.e., the tenth column of the matrix, including the column of rsids).  The second column is the "pattern type," which indicates the allele combination when a single fixed chromosome is placed next to the fixed one.  The third column is the list of SNPs (starting from 01) at which the indiciated pattern is observed.  Chromosome Overlap will produce a <i>meta-chromosome</i> <kbd>01_0,02_0,04_0,07_0,09_1</kbd> indicating SNPs 1, 2, 4, and 5 are all in the REF state and SNP 9 is in the ALT state; this is the pattern shared by the fixed chromosome and the variable chromosome.</p>

<p>All meta-chromosomes are collected with their counts in the file <kbd>Pattern.NAME.2,j.txt</kbd>. The name <kbd>2,j</kbd> indicates a fixed chromosome <kbd>2</kbd> has been combined with all variable chromosomes <kbd>j</kbd>.  But since we want to range through all fixed chromosomes and save the patterns in one place, the patterns generated by <kbd>3,j</kbd>, <kbd>4,j</kbd>, etc. have been collected here as well.  If you had indicated that a range of tuples of fixed chromosomes (i.e., columns of the transposed matrix) be considered, say 2 to 51 and 52 to 101, another file <kbd>Pattern.NAME.52,j.txt</kbd> would have been generated.  You would then combine different Pattern files using</p>

<pre>
<code>
sh HOME_DIR/ChromosomeOverlap_initiation_combine.sh NAME 1 "Iteration000"
</code>
</pre>

<p><kbd>NAME</kbd> is the name variable you used in the previous step to name your output file.  We will not need to deal with the second input, but it should be left as "1."  The name Iteration000 will be appended to the combined file so that we know these meta-chromosomes have come from the initial round.  This program workds by finding the "minimum" label by computing the "distance" between each number (e.g., the 2 in <kbd>2,j</kbd> and the 50 in <kbd>50,j</kbd>) and its position in the string (1 for both); clearly this favors suffixes with small numbers in the leftmost positions.  The combined file <kbd>Pattern_combined.Iteration000.NAME.2,j.txt</kbd> is appended with the minimum index over all Pattern files which have been combined.  The files used to make the combined file will be deleted.</p>

<p>A note is needed on the other files you may see.  For any combination of <var>&sigma;</var> chromosomes <var>i<sub>1</sub>,...,i<sub>&sigma;-1</sub></var> there will be several different paritions of the chromosomes into non-empty subsets.  The number of paritions of <var>&sigma;</var> things into <var>k</var> non-empty subsets is given by the <i>Stirling number</i> <var>S<sub>&sigma;,k</sub></var> (see the function <kbd>StirlingSum.sh</kbd>).  Chromosome Overlap currently only considers the single partition into one subset.  But for two chromosomes it is also possible to divide the combination into two subsets of one chromosome each; the alleles that are exclusive to each chromosome are output as a <i>split pattern</i> in the file appended <kbd>2+j</kbd>; you'll notice that these patterns are complementary.  You may notice some intermediate files with names like <kbd>2,j_2,j_2,j</kbd> or <kbd>2,j_2+j_j</kbd>, which indicate (i) the fixed/variable chromosomes being considered (i.e., <var>2</var> and <var>j</var>), (ii) the partition being evaluated (<var>2,j</var> for one subset or <var>2+j</var> for two subsets), and (iii) the meta-chromosomes exclusive to that subset (e.g., the alleles shared by 2 and <var>j</var> or by <var>j</var> alone.  Larger values of <var>&sigma;</var> can be used, but they will produce more subsets. The second input in <kbd>ChromosomeOverlap_initiation_combine.sh</kbd> indicates whether (1) or not (0) partitions into the same number of subsets should be combined into one output file.</p>  
 
 <h2>Filtering</h2>
 
 <p>It is not practical to iteratively overlap all the meta-chromosomes from the initiation step.  The next steps run Fisher's exact test on all the meta-chromosomes and determine those to keep based on a threshold.  You will need haplotypes files for two populations, e.g., cases and controls, with the same snp columns, as shown above.  The first script</p>
 
 <pre>
 <code>
 sh ChromosomeOverlap_haplotype_count.v2.sh cases_haplotypes,controls_haplotypes Pattern_combined.Iteration000.NAME.2,j.txt "" Iteration000.NAME
 </code>
 </pre>
 
 <p>creates a file called <kbd>haplotype_segregation.Iteration000.NAME.patterns_0000-XXXX.txt</kbd> with the meta-chromosome in the first column, the number <var>a</var> of case chromosomes with the pattern, the number <var>b</var> of case chromosomes without the pattern, the number <var>c</var> of controls chromosomes with the pattern, and the number <var>d</var> of controls chromosomes without the pattern.  The skipped input is for specifying a range of haplotypes to evaluate, so that the work can be split up across different cores (see below).</p>
 
 <p>The data in this file are sufficient to run Fisher's exact test in R:</p>
 
 <pre>
 <code>
 Rscript ChromosomeOverlap_fisher_exact.R -f haplotype_segregation.Iteration000.NAME.patterns_0000-XXXX.txt -o NAME
 </code>
 </pre>
 
 <p>creates a file <kbd>fisher_exact.Iteration000.NAME.patterns_0000-XXXX.txt</kbd> in your directory with pattern name, cases frequency, controls frequency, odds ratio, and p-value for all XXXX+1 patterns.</p>
 
 <p>It is then a simple matter to filter patterns in <kbd>Pattern_combined.Iteration000.NAME.2,j.txt</kbd> for use in the overlap iterations.  A good strategy is to rename the <kbd>Pattern_combined.Iteration000.NAME.2,j.txt</kbd> file <kbd>Pattern_combined_<b>old</b>.Iteration000.NAME.2,j.txt</kbd> with</p>
 
 <pre>
 <code>
 mv Pattern_combined.Iteration000.NAME.2,j.txt Pattern_combined_old.Iteration000.NAME.2,j.txt
 </code>
 </pre>
 
 and then run</p>
 
 <pre>
 <code>
 awk &lsquo;(NR==FNR && $5<0.05){hap[$1]; next} ($1 in hap){print $0}&rsquo; fisher_exact.NAME.txt Pattern_combined_old.Iteration000.NAME.2,j.txt > Pattern_combined.Iteration000.NAME.2,j.txt
 </code>
 </pre>
 
 <p>to create a shorter list of patterns for the iterations, each with <var>p < 0.05</var>.  It is recommended that you start with fewer than 100 patterns (and sometimes even less), or else the overlaps will quickly blow up (i.e., generate more than a few hundred thousand unique patterns that each have to be overlapped with every other).  As long as you have the original list <kbd>Pattern_combined_old.Iteration000.NAME.2,j.txt</kbd> of patterns, you can experiment with different filtering thresholds.</p>
 
 <h2>Iteration</h2>

<p>Once you have the list of filtered patterns, it's time to iterate.  The idea is to iteratively form all combinations that yeild a non-null overlap.  This could in theory be accomplished by looping through all <var>2<sup>XXXX+1</sup></var> patterns, but you would be waiting for a long time.  Instead we form all pairs of patterns at the beginning of Iteration001 and look for unique patterns.  Then we form all pairs of these patterns at the beginning of Iteration002.  The process continues until there are no more patterns, because no non-null overlaps can be generated.  This process can still blow up, but in an HPC environment it can be completed for reasonable-length patterns (~100 SNPs) and a modest number of starting patterns.  Alternatively, it is possible to terminate the process early when most of the unique patterns have already appeared.</p>

<p>To start the process, run</p>

<pre>
<code>
sh HOME_DIR/ChromosomeOverlap_iteration_sub.sh NAME 2 2,j "" "DIRECTORY" HOME_DIR
</code>
</pre>

<p>If <kbd>HOME_DIR</kbd> is on your <kbd>PATH</kbd>, you do not need to enter the entire line.  <kbd>NAME</kbd> is the identifier we have been using; it should direclty follow "Iteration000" in your <kbd>Pattern_combined</kbd> file.  The next entries indicate <var>&sigma;=2</var> and the pattern type <var>2,j</var> where a fixed chromosome was overlapped with a variable chromosome and the intersection taken.  You will usually not need to change these values.  The missing entry is the number of tuples computed per job, which is automatically altered to maintain a fixed number (100) of jobs.  Since we will be running the jobs in serial, we will leave this parameter alone.</p>

<p>The overlaps are actually evaluated using the function</p>

<pre>
<code>
sh HOME_DIR/ChromosomeOverlap_iteration.sh Pattern_combined.Iteration00X.NAME.txt "" 1000000 2 "Iteration00(X+1)" "DIRECTORY" "HOME_DIR"
</code>
</pre>

<p>This function overlaps a range of the <var>&sigma;</var>-tuples of rows of the input file (actually columns of the transposed input file).  The second input is automatically filled in and determines what range of tuples to do using the <kbd>index2combo2.sh</kbd> function.  The third input is fixed an simply means that at most 1000000 tuples will be evaluated before another step will be started.  As before, <var>&sigma;=2</var>, and <kbd>Iteration00(X+1)</kbd> names the next series of overlaps.  You will notice that this function creates a temporary sub-directory Iteration00X.Step01 in which a series of files called <kbd>Overlap_tuples....</kbd> will be generated.  These get combined by the function <kbd>ChromosomeOverlap_iteration_combine</kbd> into the new <kbd>Pattern_combined.Iteration00(X+1)</kbd> file and are then deleted.</p>

<p>Another function of the submission script is to find patterns that <i>appear</i> or <i>disappear</i> at the end of each iteration.  Another script <kbd>file_compare.sh</kbd> simply looks for lines in <kbd>TEST_FILE</kbd> that are not in <kbd>REF_FILE</kbd>.  When <kbd>TEST_FILE</kbd> is <kbd>Pattern_combined.Iteration00(X+1)</kbd> and <kbd>REF_FILE</kbd> is <kbd>Pattern_combined.Iteration00X</kbd>, the new patterns that have appeared at the end of iteration X+1 are produced (i.e., <i>pre-closed patterns</i>); when <kbd>TEST_FILE</kbd> is <kbd>Pattern_combined.Iteration00X</kbd> and <kbd>REF_FILE</kbd> is <kbd>Pattern_combined.Iteration00(X+1)</kbd>, the old patterns that have disappeared at the end of iteration X+1 are produced (i.e., <i>post-closed patterns</i>).  Because every pattern that appears must subsequently disappear, and the pre-closed patterns should eventually become post-closed patterns.  The pre- and post-patterns are tabulated, together with the iteration of the appearance or disappearance, in <kbd>Closed_patterns_pre.NAME.stats</kbd> and <kbd>Closed_patterns_post.NAME.stats</kbd>, respectively.  These closed patterns may be used in downstream analysis.</p>

<h2>A note on parallelization</h2>

<p>When an HPC cluster is available, the <kbd>_sub</kbd> jobs can be modified to submit the jobs <kbd>ChromosomeOverlap_initiaion.sh</kbd>, <kbd>ChromosomeOverlap_haplotype_count.sh</kbd>, and <kbd>ChromosomeOverlap_interation.sh</kbd> in parallel.  LSF allows users to submit an array of jobs indexed by a parameter called <kbd>$LSB_JOBINDEX</kbd>.  Each of the three functions takes a parameter <kbd>RANGE</kbd> which is of the form <kbd>STEP_SIZE.STEP_NO</kbd>.  By picking a suitable <kbd>STEP_SIZE</kbd> (determined by the number of jobs you can run at once), each job evaluates a different <kbd>STEP_NO</kbd>.  The outputs of the different jobs are combined by the functions <kbd>ChromosomeOverlap_initiaion_combine.sh</kbd> and <kbd>ChromosomeOverlap_interation_combine.sh</kbd>, so the only additional programming needed is the initiation of the job array.</p>

<p>Three functions are supplied for users with access to LSF.  To initiate the overlap process, run (or submit using <kbd>bsub</kbd>)</p>

<pre>
<code>
sh HOME_DIR/ChromosomeOverlap_initiation_bsub.sh transpose_file 1 "NAME" "50" "DIRECTORY" "HOME_DIR"
</code>
</pre>

<p>Here the fourth input "50" is the initial step size, being the number of overlaps to run in one job.  Depending on the number of tuples of fixed chromosomes, the initial step size will be multipled by 2 until the total number of jobs does not exceed 100.  This submission script will run <kbd>ChromosomeOverlap_initiation_sub.sh</kbd>, <kbd>ChromosomeOverlap_initiation.sh</kbd>, and <kbd>ChromosomeOverlap_initiation_combine.sh</kbd>.</p>

<p>After the initiation step, Fisher's exact test <i>p</i>-values can be computed by running or submitting</p>

<pre>
<code>
sh HOME_DIR/ChromosomeOverlap_haplotype_count_sub.v2.sh cases_haplotypes,controls_haplotypes Pattern_combined.Iteration000.NAME.2,j.txt 50 "" "Iteration000.NAME" DIRECTORY HOME_DIR
</code>
</pre>

<p>This script gets counts of the haplotype patterns in <kbd>Pattern_combined.Iteration000.NAME.2,j.txt</kbd> in each of <kbd>cases_haplotypes</kbd> and <kbd>controls_haplotypes</kbd> using <kbd>ChromosomeOverlap_haplotype_count.R</kbd> in groups of 50 (or more) haplotypes per job (such that the total number of jobs does not exceed 100).  The raw counts are then fed into <kbd>ChromosomeOverlap_fisher_exact.R</kbd> to get the <i>p</i>-values, which are then combined into a single file called <kbd>fisher_exact.Iteration000.NAME.patterns_0001-XXXX.txt</kbd> containing <var>XXXX</var> to be used for filtering.  The skipped input is to specify an iteration from which to select patterns for testing; the files <kbd>Closed_patterns_pre.NAME.stats</kbd> and <kbd>Closed_patterns_post.NAME.stats</kbd>, for example, which you generate during the iteration steps have iteration number in their first column to be filtered on.  The bsub script loads and R module in order to run the <kbd>Rscript</kbd> command, a process which may not work on all systems; in this case you will need to run the command manually on the <kbd>haplotype_segregation.NAME.pattern_YYYY-XXXX.txt</kbd> files (see above) and combine the results.</p>

<p>After filtering, the iteration step is almost exactly the same as before.  Run or submit</p>

<pre>
<code>
sh HOME_DIR/ChromosomeOverlap_iteration_sub_parallel.sh NAME 2 2,j "50" "DIRECTORY" HOME_DIR
</code>
</pre>

<p>This script performs the same steps as <kbd>ChromosomeOverlap_iteration_sub.sh</kbd>, except that the calls to <kbd>ChromosomeOverlap_iteration.sh</kbd> are submitted as a job array.  Again, the number over <var>&sigma;</var>-tuples of meta-chromosomes submitted per jobs is initiallized at 50 and incremented until the total number of jobs does not exceed 100.  With parallelization it is possible to submit a file of several hundred thousand meta-chromosomes on ~100-200 total SNPs.</p>

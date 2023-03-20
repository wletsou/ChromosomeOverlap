# Chromosome Overlap #

This document outlines how to perform Chromosome Overlap in an HPC environment (LSF).  The goal is to obtain a list of closed haplotype patterns from (filtered meta-)chromosomes from a population.

### Initiation ###

Your input should be a haplotype file (<kbd>.txt</kbd>) of the following form:

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
 
Each row is chromosome from subject sid; note that each has two chromosomes with the same id.  Alternate alleles of SNP </kbd>rsid_Alt</kbd> are indicated with a 1 and reference alleles with a 0.

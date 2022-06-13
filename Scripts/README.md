*** Instructions ***

Your input should be a haplotype file of the following form:

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
 
 Each row is chromosome from subject sid; note that each has two chromosomes with the same id.  Alternate alleles of SNP rsid_Alt are indicated with a 1 and reference alleles with a 0.  Transpose the file with one of the following commands.  Note that the sid row is purposely excluded.
 
 <pre>
 <code>
 n_fields=$(awk 'END{print NF}' file)
 for i in `seq 2 $n_fields`
 do
    cut -f$i file | paste -s # paste column i as a row
 done > transpose_file
 </code>
 </pre>
 
    
    

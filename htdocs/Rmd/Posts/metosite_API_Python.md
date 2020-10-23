
# How to access the API of MetOSite using Python

Python is one of the most widely used  programming languages by data scientiest. Therefore, Python is a fine choice to analyse the data supplied by *MetOSite*. In this tutorial we will show how to pull data into the Python environment from our database (DB) using its API.

The *MetOSite*'s API offers a number of end-point functions. However, the procedure to invoke them and parse the data they will return is similar in all of them:

* Make a "GET" request to pull raw data into our Python environment.
* Parse raw data through JavaScript Object Notification (JSON) into a usable format.

## Getting a summary of MetOSite data

To illustrate this two-steps procedure with a straightforward example, let us suppose we want to obtain a dataframe containing a summary (statistics) of the data found into *MetOSite*: what species are represented and with how many proteins and with how many MetO sites contributes each of these species.

For this task we will use the *Request* library:


```python
import requests
```

Afterwards, we are ready to make our first "GET" request:


```python
if  __name__ == '__main__':
    
    call = 'https://metosite.uma.es/api/summaries/species' # The API url
    response = requests.get(call)
    
    if response.status_code == 200:   
        json_species = response.json() # This is a list of dictionary
    else:
        print(response)
```

The object we have called *json_species* is a list of dictionaries, where each dictionary has three keys (*species*, *proteins* and *sites*) and three values (the species' name, the number of proteins and the number of sites). This list can be easily converted to a dataframe using the *Pandas* library:


```python
import pandas as pd
df_species = pd.DataFrame(json_species)
```

So, this *df_species* we have created is just the dataframe that we were looking for from the beginning.


```python
df_species.set_index('species').head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>proteins</th>
      <th>sites</th>
    </tr>
    <tr>
      <th>species</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Arabidopsis thaliana</th>
      <td>386</td>
      <td>530</td>
    </tr>
    <tr>
      <th>Aspergillus nidulans</th>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Bacillus cereus</th>
      <td>219</td>
      <td>500</td>
    </tr>
    <tr>
      <th>Bacillus licheniformis</th>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Bos taurus</th>
      <td>2</td>
      <td>10</td>
    </tr>
  </tbody>
</table>
</div>



Alternatively, we can get a summary putting the focus on the oxidants instead of the species. In this case the end-point function from the API that we have to call is not *species* but *oxidant*:


```python
if  __name__ == '__main__':
    
    call = 'https://metosite.uma.es/api/summaries/oxidant' # The API url
    response = requests.get(call)
    
    if response.status_code == 200:   
        json_oxidants = response.json() # This is a list of dictionary
    else:
        print(response)

df_oxidants = pd.DataFrame(json_oxidants)
df_oxidants.set_index('oxidant').head(6)
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sites</th>
    </tr>
    <tr>
      <th>oxidant</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Lipid hydroperoxides (LOOH)</th>
      <td>1</td>
    </tr>
    <tr>
      <th>NOR3</th>
      <td>1</td>
    </tr>
    <tr>
      <th>chloramine-T</th>
      <td>7</td>
    </tr>
    <tr>
      <th>peroxinitrite (ONOO)</th>
      <td>1</td>
    </tr>
    <tr>
      <th>high light irradiation</th>
      <td>528</td>
    </tr>
    <tr>
      <th>hydrogen peroxide (H2O2)</th>
      <td>2572</td>
    </tr>
  </tbody>
</table>
</div>



## Getting all the sites involved in PPI effects

Now that we have gained confidence in our abilities to access *MetOSite* through its API, we can face a task a bit more elaborated.

Basically, what we want to do is to filter the DB to keep only those entries related to changes in a biological property such as the ability to stabilize/destabilize protein-protein interactions. In order to understand how the API filters the DB using the end-point *search*, we have to introduce previously some basic ideas related to the *Groups* and *Functional Categories*.

Each MetO site is assigned to one of three possible *Groups*. **Group 1** is composed of all those MetO sites coming from high-throughput studies, reason for which nothing is known about the effect of their sulfoxidation on the protein's properties, just because it has not been addressed. On the contrary, the effect of the oxidation of residues belonging to **Group 2** has been assessed, but no effect could be found. Finally, **Group 3** encompasses all the methionine sites whose sulfoxidation has been reported to have an effect on at least one of the six following biological properties:

* Gain of activity
* Loss of activity
* Gain of protein-protein interaction
* Loss of protein-protein interaction
* Effect on protein stability
* Effect on subcellular localization

Each of these six properties can be considered as a binary variable. Thus, a value of 1 for any of them means that experimental evidence supporting such an effect has been published. On the contrary, a value of 0 only means that we have not found experimental evidence to support such an effect. In this way, attending to these variables, we have <img src="https://latex.codecogs.com/svg.latex?2^6" /> *Functional Categories* (FCs). In other words, the FC of a given MetO site can be encoded by a vector of dimension 6. It should be noted that a site with a vector (0,0,0,0,0,0), meaning that no effect has been described for the oxidation of that site, can belong to *Group 1* or *Group 2*. That is, actually we will deal with <img src="https://latex.codecogs.com/svg.latex?2^6+1=65" /> FCs.

### Making use of the *mapping*  end-point function

Thus, to know what FCs correspond with sites involved in changes of protein-protein interaction (PPI), we will make use of an ancillary end-point function called *mapping*. This function takes two arguments. The first one is related to the *Group*. For instance, the string *001* is interpreted as we are only interested in *Group 3*. If instead of *001* we pass *101* as the first argument, then *mapping* will interprets that we want to filter out *Group 2* and kepp the groups 1 and 3. The second argument that should be passed to *mapping* is a 6-dimensional vector that provides information about the effect on the six biological properties listed above. For instance, with the point (0,0,1,1,0,0) we would retrieve those sites for which a gain and loss of PPI has been reported (probably with different partners), but no other effect on the remaining properties has been described. On the other hand, if we have a site causing a gain of PPI and a loss of PPI but also a gain of activity, then the rigth argument would be (1,0,1,1,0,0). Please, note that the first and second coordinates point to the gain and loss of activity, respectively, and so on (keeping the order of the list given above).

What if we are interested in those sites that when oxidized lead to a gain and a loss of PPI but we do not care about the four remaining properties (that may or may not be affected). In this case, the right argument would be (2,2,1,1,2,2). As the insightful reader would have intuited, the integer 2 means: it does not matter whether the property has been described to be affected or not. For instance, (2,0,1,1,0,0) = (0,0,1,1,0,0) <img src="https://latex.codecogs.com/svg.latex?\cup" /> (1,0,1,1,0,0).

In the case we are developing herein, we are interested in those MetO sites involved in gain and/or loss of PPI without any other consideration. So, these sites are encoded as follows (2,2,1,2,2,2) <img src="https://latex.codecogs.com/svg.latex?\cup" /> (2,2,2,1,2,2). At this point, we are in conditions to use knowingly the *mapping* end-point


```python
groups = '001' # only Group 3 is goint to be relevant for us

## ------------------------- Sites gaining PPI ----------------------------- ##

categories = '221222' # Note we use neither parenthesis nor commas.
call = 'https://metosite.uma.es/api/sites/mapping/'+groups+'/'+categories

if  __name__ == '__main__':
    
    response = requests.get(call)
    
    if response.status_code == 200:   
        gPPI = response.json() # list containing the requested FCs
    else:
        print(response)

## -------------------------- Sites losing PPI ----------------------------- ##

categories = '222122' # Note we use neither parenthesis nor commas.
call = 'https://metosite.uma.es/api/sites/mapping/'+groups+'/'+categories

if  __name__ == '__main__':
    
    response = requests.get(call)
    
    if response.status_code == 200:   
        lPPI = response.json() # list containing the requested FCs
    else:
        print(response)

## ------------------------- Joining btoh sets ----------------------------- ##
fcPPI = (gPPI + lPPI)
fcPPI.sort()
fcPPI = set(fcPPI)
print(fcPPI)
print('\t')
print('There are ' + str(len(fcPPI)) + ' funtional categories we are interested in')

{5, 6, 9, 10, 11, 12, 17, 18, 19, 20, 21, 24, 25, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65}
```

    
There are 48 functional categories we are interested in.


Now that we know there are 48 different functional categories that meet the requiriments (be involved in PPI), we can move to find the MetO sites present into *MetOSite* that belong to any of these FCs.

### Making use of the *search* end-point function

This function takes three arguments: (i) the first one is related to the FCs we wish to retrieve,
(ii) the second one allows us to filter using a taxon criterium, that is, the organism(s) we are interested in, and (iii) the third criterium is related to the oxidant(s) we want to consider.

Because the FCs we are going to pass to the *search* function need to be separated from each other by the symbol *&*, we are going to write a function that will make the formatting work for us:


```python
def format_FC(fc):
    formatted_fc = ''
    for c in fc:
        formatted_fc = formatted_fc + str(c) + '&'   
    formatted_fc = formatted_fc[:-1]
    return(formatted_fc)
```

So, let's use that function to format the set of 48 FCs we got previously, and then make use of the *search* end-point function:


```python
ffc = format_FC(fcPPI) # formatted FCs related to PPI
organism = '-1' # meaning we don't care about the organism
oxidant = '-1' # meaning we don't care about the oxidant

call = 'https://metosite.uma.es/api/sites/search/'+ffc+'/'+organism+'/'+oxidant

if __name__ == '__main__':

    response = requests.get(call)

    if response.status_code == 200:
        json_results = response.json()
    else:
        print(response)


ppi_results = pd.DataFrame(json_results)
ppi_results.set_index('prot_name').head(5)

```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>met_pos</th>
      <th>met_vivo_vitro</th>
      <th>org_oxidant</th>
      <th>org_sp</th>
      <th>prot_id</th>
      <th>reg_id</th>
    </tr>
    <tr>
      <th>prot_name</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Amyloid beta-peptide(1-42)</th>
      <td>35</td>
      <td>vitro</td>
      <td>hydrogen peroxide (H2O2)</td>
      <td>Homo sapiens</td>
      <td>P05067</td>
      <td>10</td>
    </tr>
    <tr>
      <th>Apolipoprotein D</th>
      <td>93</td>
      <td>vitro</td>
      <td>lipid hydroperoxides (LOOH)</td>
      <td>Homo sapiens</td>
      <td>P05090</td>
      <td>12</td>
    </tr>
    <tr>
      <th>Cofilin-1</th>
      <td>115</td>
      <td>vitro</td>
      <td>taurine chloramine (TnCl)</td>
      <td>Homo sapiens</td>
      <td>P23528</td>
      <td>18</td>
    </tr>
    <tr>
      <th>alpha-Actin</th>
      <td>46</td>
      <td>both</td>
      <td>mical-catalyzed</td>
      <td>Oryctolagus cuniculus</td>
      <td>P68135</td>
      <td>18</td>
    </tr>
    <tr>
      <th>Actin-5C</th>
      <td>44</td>
      <td>both</td>
      <td>mical-catalyzed</td>
      <td>Drosophila melanogaster</td>
      <td>P10987</td>
      <td>18</td>
    </tr>
  </tbody>
</table>
</div>



Thus, the object we have named *ppi_results* is the dataframe we wanted.

## Final remarks

This tutorial does not pretend to be exhaustive, on the contrary, it aims to be a primer from which the user can continue on his/her own, exploring how to communicate with the API of *MetOSite*. We encourage the user to explore the other end-points that will find at https://metosite.uma.es/api-docs

## Further readings

For further details you may want to check the following tutorial:
https://www.digitalocean.com/community/tutorials/how-to-use-web-apis-in-python-3#step-4-%E2%80%94-working-with-a-different-api




Copyright &#169; 2018 The MetOSite team

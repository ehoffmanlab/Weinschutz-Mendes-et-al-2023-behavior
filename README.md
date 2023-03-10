# **Visual-Startle Analysis Code**
This code was used in Weinschutz Mendes et al. (2023) to analyze visual-startle phenotypes in zebrafish mutants of autism genes. This method is described in the STAR Methods section of Weinschutz Mendes et al. (2023). 
# **Required user input files**
1. ViewPoint file (.XLS)\*
1. Grouping file (.txt)\*\*
1. Timing file (.txt)\*\*

\* Generated by ZebraLab software (Viewpoint Life Sciences, Montreal, Quebec, Canada)

\*\* Generated by the user
## **ViewPoint file(s)**
- This is the .XLS file generated by ZebraLab software. 
- To analyze visual-startle responses, **ZebraLab quantization mode** (Zebrabox and ZebraLab software) is used to record larval activity every 1 second. 
- Visual-startle-off (**VSR-OFF**) and visual-startle-on (**VSR-ON**) experiments are run as separate protocols in ZebraLab. 
  - For VSR-OFF, larvae are acclimated to white light background illumination (light intensity = 15%) for 1 hour, after which their baseline activity is tracked for 30 minutes followed by five 1-second flashes of darkness (lights turn off) at 29-second intervals. 
  - For VSR-ON, the assay is reversed, where larvae were acclimated to darkness and then exposed to 1-second flashes of bright white light (100% light intensity) at 29-second intervals. 
  - The total duration of each experiment is 92.5 minutes (5,550 seconds): 60 min, acclimation; 30 min, baseline; 2.5 min, stimuli/post-stimuli (1 second stimulus, 29 seconds post-stimulus, x 5 stimuli). 
- Zebralab threshold parameters for detection are: detection threshold, 40; burst, 25; freeze, 4; bin size, 1 second. 
- To be compatible with the MATLAB analysis code, file names should use the following convention: “**YYMMDD\_0X\_off**” (VSR-OFF) and “**YYMMDD\_0X\_on**” (VSR-ON), where “xx” refers to an arbitrary Zebrabox number specified by the laboratory (e.g. 0A, 01).
## **Grouping file**
- This is a .txt file that specifies the experimental group (e.g., genotype).
- To create genotype file, include the following in a .txt file—
  - First row: Genotype 1, Genotype 2, etc…
  - Second row: names of different genotypes (or experimental groups), e.g., wild-type, heterozygous, homozygous. 
- In each column, write the fish IDs that corresponds to each genotype. 
  - The fish IDs correspond to the location on a 96-well plate, where the rows are A through H and the columns are 1 through 12: A1 = 1, A2 = 2, A3= 3, … B1 = 13, B2 = 14, B3 = 15, … H1 = 85, H2 = 86, H3 = 87 … H12 = 96.
- The genotype file must be saved in the same format as the ZebraLab file (e.g., YYMMDD\_0Xgenotype.txt) 
- Save as Tab Delineated Text (.txt)
- **Sample Genotype file syntax:**


|Genotype 1|Genotype 2|Genotype 3|Genotype 4|Genotype 5|etc...|
| :- | :- | :- | :- | :- | :- |
|HOM:HOM|HOM:HET|HOM:WT|HET:HOM|HET:HET||
|5|36|14|29|22||
|18|54|27|33|24||
|91|81|80|77|40||
|x|x|x|x|x||

## **Timing File(s)**
- This is a .txt file that specifies the timestamps of events (e.g., stimuli) in the experiment.
- For **VSR-OFF**t, the timing file is labeled as “**lightsoff\_Timing.txt**” and has the following syntax:

|Onset|Offset|Light|
| :- | :- | :- |
|1|5400|on|
|5400|5401|off|
|5401|5430|on|
|5430|5431|off|
|5431|5460|on|
|5460|5461|off|
|5461|5490|on|
|5490|5491|off|
|5491|5520|on|
|5520|5521|off|
|5521|5550|on|

For **VSR-ON**, the timing file is labeled as “**lightson\_Timing.txt**” and has the following syntax:

|Onset|Offset|Light|
| :- | :- | :- |
|1|5400|off|
|5400|5401|on|
|5401|5430|off|
|5430|5431|on|
|5431|5460|off|
|5460|5461|on|
|5461|5490|off|
|5490|5491|on|
|5491|5520|off|
|5520|5521|on|
|5521|5550|off|

## **To run the visual-startle analysis code:**
**\*\*\* MATLAB 2019 is required to run this code.**

1. Add matlab files to path: The first time the program is run on a computer, the location of the .m files must be added to the path.
1. The user will be asked to specify the name of the project, the save location for output files, and the locations of the grouping, Viewpoint, and timing files. See below--

projectName = 'NAME';

saveLocation = '/Users/...';



groupingFile = {'/Users/...'};



viewPointFile = {'/Users/...'};



timingFile = '/Users/...';

3. Run “**ABA\_CTRL.m” in the “animalBehavioranalysis” folder.**
## **Outputs**
1. Data object .m file
   1. The Data object is a complex data type called an “object.” 
   1. The Data Object contains the raw experiment data, metadata information (such as grouping specifications), specific processed results, and methods for manipulating the data.  As such, it is the main tool researchers will use to process and report data.
   1. After constructing the internal data object, the program saves a copy of this object and a Matlab .m file (format is **YYMMDD\_0x\_off.mat** for VSR-OFF or **YYMMDD\_0x\_on.mat** for VSR-ON). This file can be loaded at a later point and contains both the raw data, processed results, and methods for analyzing the data. 
1. Results folder
   1. Aggregated raw data file
      1. After reading relevant input files, the program builds an organized internal data object that integrates and reconciles the disparate raw inputs. 
      1. The program writes this aggregated data to a file so that researchers have the option of examining the raw data in their preferred statistical analysis software. 
   1. Time series graphs
   1. One-way ANOVA analyses and graphs

[![DOI](https://zenodo.org/badge/601389493.svg)](https://zenodo.org/badge/latestdoi/601389493)

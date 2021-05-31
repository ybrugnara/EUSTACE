This software was developed by Antonello Squintu (KNMI) in the framework of the EUSTACE project (https://www.eustaceproject.org/).

To run the software:

1. Install required R packages (gdata and trend)
2. Prepare input files (follow instructions in the header of file AA_quant_match_MO.R or use script format_input.R)
3. Edit the input parameters in AA_quant_match_MO.R (if needed)

If only one series has to be homogenized:
4. Uncomment and edit the paths in AA_quant_match_MO.R
5. Run AA_quant_match_MO.R

More than one series:
4. Edit the input parameters in run.R
5. Run run.R
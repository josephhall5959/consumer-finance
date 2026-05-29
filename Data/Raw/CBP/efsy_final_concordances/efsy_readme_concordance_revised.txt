Using the industry concordances:

Each concordance has four columns:
(input ind code) -- Industry code vintage to be converted
(output ind code)-- Corresponding industry code in subsequent vintage
mappings		-- Count of number of mappings
wt_mappings	-- Mapping weight

Download the concordance corresponding to the correct CBP year.
    Years	|   Input industry code
----------------+----------------------------
1977-1987	| SIC77
1988-1997	| SIC87
1998-2002	| NAICS97
2003-2007	| NAICS02
2008-2011	| NAICS07

To convert between successive industry code vintages, merge the concordance to the data 
by the (orig ind code), multiply the original employment count by the weight to 
get the (output ind code) employment, and collapse by FIPS code and (output ind code). Concordances can be joined together to convert to later vintages.

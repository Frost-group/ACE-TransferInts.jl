begin
using Pkg
Pkg.activate(".")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("DelimitedFiles")
Pkg.add("LaTeXStrings")
Pkg.add("MultivariateStats")
Pkg.add("Plots")
Pkg.add("Suppressor")
Pkg.add("ExtXYZ")
Pkg.add("JuLIP")
Pkg.add("Plots")
Pkg.add("LaTeXStrings")

using CSV
using DataFrames
using DelimitedFiles
using JuLIP
using Plots
using LaTeXStrings
using Plots
using LaTeXStrings
end


function read_J(namefile)
	data_J = (CSV.File(namefile, delim=" ")|> Tables.matrix)
	para, J, logJ = data_J[:, 1:6], data_J[:, 7], data_J[:, 8]
	return para, J, logJ
end

para, J², logJ² = read_J("out.dat");
para[1:49,:]

J = zeros((49*97))
begin
#____1
j = 0
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 20 + j * 49] = +sqrt.(J²[05 + j * 49: 20 + j * 49])
J[21 + j * 49: 29 + j * 49] = -sqrt.(J²[21 + j * 49: 29 + j * 49])
J[30 + j * 49: 45 + j * 49] = +sqrt.(J²[30 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____2
j = 1
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 20 + j * 49] = +sqrt.(J²[05 + j * 49: 20 + j * 49])
J[21 + j * 49: 28 + j * 49] = -sqrt.(J²[21 + j * 49: 28 + j * 49])
J[29 + j * 49: 45 + j * 49] = +sqrt.(J²[29 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____3
j = 2
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 20 + j * 49] = +sqrt.(J²[05 + j * 49: 20 + j * 49])
J[19 + j * 49: 28 + j * 49] = -sqrt.(J²[19 + j * 49: 28 + j * 49])
J[29 + j * 49: 45 + j * 49] = +sqrt.(J²[29 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____4
j = 3
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 11 + j * 49] = +sqrt.(J²[05 + j * 49: 11 + j * 49])
J[12 + j * 49: 27 + j * 49] = -sqrt.(J²[12 + j * 49: 27 + j * 49])
J[28 + j * 49: 45 + j * 49] = +sqrt.(J²[28 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____5
j = 4
J[01 + j * 49: 02 + j * 49] = -sqrt.(J²[01 + j * 49: 02 + j * 49])
J[03 + j * 49: 07 + j * 49] = +sqrt.(J²[03 + j * 49: 07 + j * 49])
J[08 + j * 49: 21 + j * 49] = -sqrt.(J²[08 + j * 49: 21 + j * 49])
J[22 + j * 49: 27 + j * 49] = +sqrt.(J²[22 + j * 49: 27 + j * 49])
J[28 + j * 49: 29 + j * 49] = -sqrt.(J²[28 + j * 49: 29 + j * 49])
J[30 + j * 49: 45 + j * 49] = +sqrt.(J²[30 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____6
j = 5
J[01 + j * 49: 04 + j * 49] = +sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 19 + j * 49] = -sqrt.(J²[05 + j * 49: 19 + j * 49])
J[20 + j * 49: 26 + j * 49] = +sqrt.(J²[20 + j * 49: 26 + j * 49])
J[27 + j * 49: 35 + j * 49] = -sqrt.(J²[27 + j * 49: 35 + j * 49])
J[36 + j * 49: 44 + j * 49] = +sqrt.(J²[36 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = -sqrt.(J²[45 + j * 49: 49 + j * 49])
#____7
j = 6
J[01 + j * 49: 04 + j * 49] = +sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 18 + j * 49] = -sqrt.(J²[05 + j * 49: 18 + j * 49])
J[19 + j * 49: 26 + j * 49] = +sqrt.(J²[19 + j * 49: 26 + j * 49])
J[27 + j * 49: 41 + j * 49] = -sqrt.(J²[27 + j * 49: 41 + j * 49])
J[42 + j * 49: 43 + j * 49] = +sqrt.(J²[42 + j * 49: 43 + j * 49])
J[44 + j * 49: 49 + j * 49] = -sqrt.(J²[44 + j * 49: 49 + j * 49])
#____8
j = 7
J[01 + j * 49: 04 + j * 49] = +sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 18 + j * 49] = -sqrt.(J²[05 + j * 49: 18 + j * 49])
J[18 + j * 49: 26 + j * 49] = +sqrt.(J²[18 + j * 49: 26 + j * 49])
J[27 + j * 49: 45 + j * 49] = -sqrt.(J²[27 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
#____9
j = 8
J[01 + j * 49: 04 + j * 49] = +sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 15 + j * 49] = -sqrt.(J²[05 + j * 49: 15 + j * 49])
J[16 + j * 49: 25 + j * 49] = +sqrt.(J²[16 + j * 49: 25 + j * 49])
J[26 + j * 49: 44 + j * 49] = -sqrt.(J²[26 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = +sqrt.(J²[45 + j * 49: 49 + j * 49])
#____10
j = 9
J[01 + j * 49: 04 + j * 49] = +sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 14 + j * 49] = -sqrt.(J²[05 + j * 49: 14 + j * 49])
J[15 + j * 49: 25 + j * 49] = +sqrt.(J²[15 + j * 49: 25 + j * 49])
J[26 + j * 49: 44 + j * 49] = -sqrt.(J²[26 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = +sqrt.(J²[45 + j * 49: 49 + j * 49])
end;

begin
#____11
j = 10
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 12 + j * 49] = -sqrt.(J²[06 + j * 49: 12 + j * 49])
J[13 + j * 49: 25 + j * 49] = +sqrt.(J²[13 + j * 49: 25 + j * 49])
J[26 + j * 49: 44 + j * 49] = -sqrt.(J²[26 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = +sqrt.(J²[45 + j * 49: 49 + j * 49])
#____12
j = 11
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 12 + j * 49] = +sqrt.(J²[06 + j * 49: 12 + j * 49])
J[13 + j * 49: 25 + j * 49] = +sqrt.(J²[13 + j * 49: 25 + j * 49])
J[26 + j * 49: 44 + j * 49] = -sqrt.(J²[26 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = +sqrt.(J²[45 + j * 49: 49 + j * 49])
#____13
j = 12
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 12 + j * 49] = +sqrt.(J²[05 + j * 49: 12 + j * 49])
J[13 + j * 49: 25 + j * 49] = +sqrt.(J²[13 + j * 49: 25 + j * 49])
J[26 + j * 49: 43 + j * 49] = -sqrt.(J²[26 + j * 49: 43 + j * 49])
J[44 + j * 49: 49 + j * 49] = +sqrt.(J²[44 + j * 49: 49 + j * 49])
#____14
j = 13
J[01 + j * 49: 05 + j * 49] = -sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 12 + j * 49] = +sqrt.(J²[06 + j * 49: 12 + j * 49])
J[13 + j * 49: 24 + j * 49] = +sqrt.(J²[13 + j * 49: 24 + j * 49])
J[25 + j * 49: 43 + j * 49] = -sqrt.(J²[25 + j * 49: 43 + j * 49])
J[44 + j * 49: 49 + j * 49] = +sqrt.(J²[44 + j * 49: 49 + j * 49])
#____15
j = 14
J[01 + j * 49: 06 + j * 49] = -sqrt.(J²[01 + j * 49: 06 + j * 49])
J[07 + j * 49: 12 + j * 49] = +sqrt.(J²[07 + j * 49: 12 + j * 49])
J[13 + j * 49: 24 + j * 49] = +sqrt.(J²[13 + j * 49: 24 + j * 49])
J[25 + j * 49: 42 + j * 49] = -sqrt.(J²[25 + j * 49: 42 + j * 49])
J[43 + j * 49: 49 + j * 49] = +sqrt.(J²[43 + j * 49: 49 + j * 49])
#____16
j = 15
J[01 + j * 49: 07 + j * 49] = -sqrt.(J²[01 + j * 49: 07 + j * 49])
J[08 + j * 49: 12 + j * 49] = +sqrt.(J²[08 + j * 49: 12 + j * 49])
J[13 + j * 49: 23 + j * 49] = +sqrt.(J²[13 + j * 49: 23 + j * 49])
J[24 + j * 49: 41 + j * 49] = -sqrt.(J²[24 + j * 49: 41 + j * 49])
J[42 + j * 49: 49 + j * 49] = +sqrt.(J²[42 + j * 49: 49 + j * 49])
#____17
j = 16
J[01 + j * 49: 08 + j * 49] = -sqrt.(J²[01 + j * 49: 08 + j * 49])
J[09 + j * 49: 12 + j * 49] = +sqrt.(J²[09 + j * 49: 12 + j * 49])
J[13 + j * 49: 23 + j * 49] = +sqrt.(J²[13 + j * 49: 23 + j * 49])
J[23 + j * 49: 41 + j * 49] = -sqrt.(J²[23 + j * 49: 41 + j * 49])
J[42 + j * 49: 49 + j * 49] = +sqrt.(J²[42 + j * 49: 49 + j * 49])
#____18
j = 17
J[01 + j * 49: 09 + j * 49] = -sqrt.(J²[01 + j * 49: 09 + j * 49])
J[10 + j * 49: 12 + j * 49] = +sqrt.(J²[10 + j * 49: 12 + j * 49])
J[13 + j * 49: 23 + j * 49] = +sqrt.(J²[13 + j * 49: 23 + j * 49])
J[23 + j * 49: 40 + j * 49] = -sqrt.(J²[23 + j * 49: 40 + j * 49])
J[41 + j * 49: 49 + j * 49] = +sqrt.(J²[41 + j * 49: 49 + j * 49])
#____19
j = 18
J[01 + j * 49: 11 + j * 49] = -sqrt.(J²[01 + j * 49: 11 + j * 49])
J[12 + j * 49: 23 + j * 49] = +sqrt.(J²[12 + j * 49: 23 + j * 49])
J[22 + j * 49: 38 + j * 49] = -sqrt.(J²[22 + j * 49: 38 + j * 49])
J[39 + j * 49: 49 + j * 49] = +sqrt.(J²[39 + j * 49: 49 + j * 49])
#____20
j = 19
J[01 + j * 49: 15 + j * 49] = -sqrt.(J²[01 + j * 49: 15 + j * 49])
J[16 + j * 49: 21 + j * 49] = +sqrt.(J²[16 + j * 49: 21 + j * 49])
J[21 + j * 49: 35 + j * 49] = -sqrt.(J²[21 + j * 49: 35 + j * 49])
J[36 + j * 49: 49 + j * 49] = +sqrt.(J²[36 + j * 49: 49 + j * 49])
end

begin
#____21
j = 20
J[01 + j * 49: 20 + j * 49] = -sqrt.(J²[01 + j * 49: 20 + j * 49])
J[21 + j * 49: 29 + j * 49] = -sqrt.(J²[21 + j * 49: 29 + j * 49])
J[30 + j * 49: 49 + j * 49] = +sqrt.(J²[30 + j * 49: 49 + j * 49])
#____22
j = 21
J[01 + j * 49: 19 + j * 49] = -sqrt.(J²[01 + j * 49: 19 + j * 49])
J[20 + j * 49: 49 + j * 49] = +sqrt.(J²[20 + j * 49: 49 + j * 49])
#____23
j = 22
J[01 + j * 49: 15 + j * 49] = -sqrt.(J²[01 + j * 49: 15 + j * 49])
J[16 + j * 49: 49 + j * 49] = +sqrt.(J²[16 + j * 49: 49 + j * 49])
#____24
j = 23
J[01 + j * 49: 49 + j * 49] = +sqrt.(J²[1 + j * 49: 49 + j * 49])
#____25
j = 24
J[1 + j * 49: 49 + j * 49] = +sqrt.(J²[1 + j * 49: 49 + j * 49])
#____26
j = 25
J[1 + j * 49: 49 + j * 49] = +sqrt.(J²[1 + j * 49: 49 + j * 49])
#____27
j = 26
J[01 + j * 49: 34 + j * 49] = +sqrt.(J²[01 + j * 49: 34 + j * 49])
J[35 + j * 49: 49 + j * 49] = -sqrt.(J²[35 + j * 49: 49 + j * 49])
#____28
j = 27
J[01 + j * 49: 30 + j * 49] = +sqrt.(J²[01 + j * 49: 30 + j * 49])
J[31 + j * 49: 49 + j * 49] = -sqrt.(J²[31 + j * 49: 49 + j * 49])
#____29
j = 28
J[01 + j * 49: 21 + j * 49] = +sqrt.(J²[01 + j * 49: 21 + j * 49])
J[22 + j * 49: 49 + j * 49] = -sqrt.(J²[22 + j * 49: 49 + j * 49])
#____30
j = 29
J[01 + j * 49: 21 + j * 49] = +sqrt.(J²[01 + j * 49: 21 + j * 49])
J[15 + j * 49: 29 + j * 49] = -sqrt.(J²[15 + j * 49: 29 + j * 49])
J[30 + j * 49: 34 + j * 49] = +sqrt.(J²[30 + j * 49: 34 + j * 49])
J[35 + j * 49: 49 + j * 49] = -sqrt.(J²[35 + j * 49: 49 + j * 49])
end

begin
#____31
j = 30
J[01 + j * 49: 11 + j * 49] = +sqrt.(J²[01 + j * 49: 11 + j * 49])
J[12 + j * 49: 28 + j * 49] = -sqrt.(J²[12 + j * 49: 28 + j * 49])
J[29 + j * 49: 39 + j * 49] = +sqrt.(J²[29 + j * 49: 39 + j * 49])
J[38 + j * 49: 49 + j * 49] = -sqrt.(J²[38 + j * 49: 49 + j * 49])
#____32
j = 31
J[01 + j * 49: 09 + j * 49] = +sqrt.(J²[01 + j * 49: 09 + j * 49])
J[10 + j * 49: 27 + j * 49] = -sqrt.(J²[10 + j * 49: 27 + j * 49])
J[28 + j * 49: 39 + j * 49] = +sqrt.(J²[28 + j * 49: 39 + j * 49])
J[40 + j * 49: 49 + j * 49] = -sqrt.(J²[40 + j * 49: 49 + j * 49])
#____33
j = 32
J[01 + j * 49: 08 + j * 49] = +sqrt.(J²[01 + j * 49: 08 + j * 49])
J[09 + j * 49: 27 + j * 49] = -sqrt.(J²[09 + j * 49: 27 + j * 49])
J[28 + j * 49: 41 + j * 49] = +sqrt.(J²[28 + j * 49: 41 + j * 49])
J[42 + j * 49: 49 + j * 49] = -sqrt.(J²[42 + j * 49: 49 + j * 49])
#____34
j = 33
J[01 + j * 49: 07 + j * 49] = +sqrt.(J²[01 + j * 49: 07 + j * 49])
J[08 + j * 49: 26 + j * 49] = -sqrt.(J²[08 + j * 49: 26 + j * 49])
J[27 + j * 49: 42 + j * 49] = +sqrt.(J²[27 + j * 49: 42 + j * 49])
J[43 + j * 49: 49 + j * 49] = -sqrt.(J²[43 + j * 49: 49 + j * 49])
#____35
j = 34
J[01 + j * 49: 07 + j * 49] = +sqrt.(J²[01 + j * 49: 07 + j * 49])
J[08 + j * 49: 26 + j * 49] = -sqrt.(J²[08 + j * 49: 26 + j * 49])
J[27 + j * 49: 43 + j * 49] = +sqrt.(J²[27 + j * 49: 43 + j * 49])
J[44 + j * 49: 49 + j * 49] = -sqrt.(J²[44 + j * 49: 49 + j * 49])
#____36
j = 35
J[01 + j * 49: 07 + j * 49] = +sqrt.(J²[01 + j * 49: 07 + j * 49])
J[08 + j * 49: 26 + j * 49] = -sqrt.(J²[08 + j * 49: 26 + j * 49])
J[27 + j * 49: 44 + j * 49] = +sqrt.(J²[27 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = -sqrt.(J²[45 + j * 49: 49 + j * 49])
#____37
j = 36
J[01 + j * 49: 06 + j * 49] = +sqrt.(J²[01 + j * 49: 06 + j * 49])
J[07 + j * 49: 25 + j * 49] = -sqrt.(J²[07 + j * 49: 25 + j * 49])
J[26 + j * 49: 45 + j * 49] = +sqrt.(J²[26 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____38
j = 37
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 25 + j * 49] = -sqrt.(J²[06 + j * 49: 25 + j * 49])
J[26 + j * 49: 49 + j * 49] = +sqrt.(J²[26 + j * 49: 49 + j * 49])
#____39
j = 38
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 24 + j * 49] = -sqrt.(J²[06 + j * 49: 24 + j * 49])
J[25 + j * 49: 38 + j * 49] = +sqrt.(J²[25 + j * 49: 38 + j * 49])
J[39 + j * 49: 44 + j * 49] = -sqrt.(J²[39 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = +sqrt.(J²[45 + j * 49: 49 + j * 49])
#____40
j = 39
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 24 + j * 49] = -sqrt.(J²[06 + j * 49: 24 + j * 49])
J[25 + j * 49: 35 + j * 49] = +sqrt.(J²[25 + j * 49: 35 + j * 49])
J[36 + j * 49: 45 + j * 49] = -sqrt.(J²[36 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
end

begin
#____41
j = 40
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 24 + j * 49] = -sqrt.(J²[06 + j * 49: 24 + j * 49])
J[25 + j * 49: 34 + j * 49] = +sqrt.(J²[25 + j * 49: 34 + j * 49])
J[35 + j * 49: 45 + j * 49] = -sqrt.(J²[35 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
#____42
j = 41
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 24 + j * 49] = -sqrt.(J²[06 + j * 49: 24 + j * 49])
J[25 + j * 49: 33 + j * 49] = +sqrt.(J²[25 + j * 49: 33 + j * 49])
J[34 + j * 49: 45 + j * 49] = -sqrt.(J²[34 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
#____43
j = 42
J[01 + j * 49: 05 + j * 49] = -sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 08 + j * 49] = +sqrt.(J²[06 + j * 49: 08 + j * 49])
J[09 + j * 49: 23 + j * 49] = -sqrt.(J²[09 + j * 49: 23 + j * 49])
J[24 + j * 49: 31 + j * 49] = +sqrt.(J²[24 + j * 49: 31 + j * 49])
J[32 + j * 49: 45 + j * 49] = -sqrt.(J²[32 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
#____44*******
j = 43
J[01 + j * 49: 05 + j * 49] = -sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 15 + j * 49] = +sqrt.(J²[06 + j * 49: 15 + j * 49])
J[16 + j * 49: 23 + j * 49] = -sqrt.(J²[16 + j * 49: 23 + j * 49])
J[24 + j * 49: 30 + j * 49] = +sqrt.(J²[24 + j * 49: 30 + j * 49])
J[31 + j * 49: 45 + j * 49] = -sqrt.(J²[31 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
#____45
j = 44
J[01 + j * 49: 05 + j * 49] = -sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 20 + j * 49] = +sqrt.(J²[06 + j * 49: 20 + j * 49])
J[21 + j * 49: 23 + j * 49] = -sqrt.(J²[21 + j * 49: 23 + j * 49])
J[24 + j * 49: 29 + j * 49] = +sqrt.(J²[24 + j * 49: 29 + j * 49])
J[30 + j * 49: 43 + j * 49] = -sqrt.(J²[30 + j * 49: 43 + j * 49])
J[44 + j * 49: 47 + j * 49] = +sqrt.(J²[44 + j * 49: 47 + j * 49])
J[48 + j * 49: 49 + j * 49] = -sqrt.(J²[48 + j * 49: 49 + j * 49])
#____46
j = 45
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 22 + j * 49] = +sqrt.(J²[05 + j * 49: 22 + j * 49])
J[23 + j * 49: 38 + j * 49] = -sqrt.(J²[23 + j * 49: 38 + j * 49])
J[39 + j * 49: 46 + j * 49] = +sqrt.(J²[39 + j * 49: 46 + j * 49])
J[47 + j * 49: 49 + j * 49] = -sqrt.(J²[47 + j * 49: 49 + j * 49])
#____47
j = 46
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 22 + j * 49] = +sqrt.(J²[05 + j * 49: 22 + j * 49])
J[23 + j * 49: 31 + j * 49] = -sqrt.(J²[23 + j * 49: 31 + j * 49])
J[32 + j * 49: 45 + j * 49] = +sqrt.(J²[32 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____48
j = 47
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 21 + j * 49] = +sqrt.(J²[05 + j * 49: 21 + j * 49])
J[22 + j * 49: 30 + j * 49] = -sqrt.(J²[22 + j * 49: 30 + j * 49])
J[31 + j * 49: 45 + j * 49] = +sqrt.(J²[31 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____49
j = 48
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 21 + j * 49] = +sqrt.(J²[05 + j * 49: 21 + j * 49])
J[22 + j * 49: 29 + j * 49] = -sqrt.(J²[22 + j * 49: 29 + j * 49])
J[30 + j * 49: 45 + j * 49] = +sqrt.(J²[30 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____50
j = 49
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 20 + j * 49] = +sqrt.(J²[05 + j * 49: 20 + j * 49])
J[21 + j * 49: 28 + j * 49] = -sqrt.(J²[21 + j * 49: 28 + j * 49])
J[29 + j * 49: 45 + j * 49] = +sqrt.(J²[29 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
end

begin
#____51
j = 50
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 18 + j * 49] = +sqrt.(J²[05 + j * 49: 18 + j * 49])
J[19 + j * 49: 28 + j * 49] = -sqrt.(J²[19 + j * 49: 28 + j * 49])
J[29 + j * 49: 45 + j * 49] = +sqrt.(J²[29 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____52
j = 51
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 11 + j * 49] = +sqrt.(J²[05 + j * 49: 11 + j * 49])
J[12 + j * 49: 27 + j * 49] = -sqrt.(J²[12 + j * 49: 27 + j * 49])
J[28 + j * 49: 45 + j * 49] = +sqrt.(J²[28 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____53
j = 52
J[01 + j * 49: 02 + j * 49] = -sqrt.(J²[01 + j * 49: 02 + j * 49])
J[03 + j * 49: 07 + j * 49] = +sqrt.(J²[03 + j * 49: 07 + j * 49])
J[08 + j * 49: 21 + j * 49] = -sqrt.(J²[08 + j * 49: 21 + j * 49])
J[22 + j * 49: 27 + j * 49] = +sqrt.(J²[22 + j * 49: 27 + j * 49])
J[28 + j * 49: 29 + j * 49] = -sqrt.(J²[28 + j * 49: 29 + j * 49])
J[30 + j * 49: 44 + j * 49] = +sqrt.(J²[30 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = -sqrt.(J²[45 + j * 49: 49 + j * 49])
#____54
j = 53
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 19 + j * 49] = -sqrt.(J²[06 + j * 49: 19 + j * 49])
J[20 + j * 49: 26 + j * 49] = +sqrt.(J²[20 + j * 49: 26 + j * 49])
J[27 + j * 49: 35 + j * 49] = -sqrt.(J²[27 + j * 49: 35 + j * 49])
J[36 + j * 49: 44 + j * 49] = +sqrt.(J²[36 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = -sqrt.(J²[45 + j * 49: 49 + j * 49])
#____55
j = 54
J[01 + j * 49: 04 + j * 49] = +sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 18 + j * 49] = -sqrt.(J²[05 + j * 49: 18 + j * 49])
J[19 + j * 49: 26 + j * 49] = +sqrt.(J²[19 + j * 49: 26 + j * 49])
J[27 + j * 49: 41 + j * 49] = -sqrt.(J²[27 + j * 49: 41 + j * 49])
J[42 + j * 49: 44 + j * 49] = +sqrt.(J²[42 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = -sqrt.(J²[45 + j * 49: 49 + j * 49])
#____56
j = 55
J[01 + j * 49: 04 + j * 49] = +sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 17 + j * 49] = -sqrt.(J²[05 + j * 49: 17 + j * 49])
J[18 + j * 49: 26 + j * 49] = +sqrt.(J²[18 + j * 49: 26 + j * 49])
J[27 + j * 49: 44 + j * 49] = -sqrt.(J²[27 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = +sqrt.(J²[45 + j * 49: 49 + j * 49])
#____57
j = 56
J[01 + j * 49: 04 + j * 49] = +sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 15 + j * 49] = -sqrt.(J²[05 + j * 49: 15 + j * 49])
J[16 + j * 49: 25 + j * 49] = +sqrt.(J²[16 + j * 49: 25 + j * 49])
J[26 + j * 49: 44 + j * 49] = -sqrt.(J²[26 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = +sqrt.(J²[45 + j * 49: 49 + j * 49])
#____58
j = 57
J[01 + j * 49: 04 + j * 49] = +sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 14 + j * 49] = -sqrt.(J²[05 + j * 49: 14 + j * 49])
J[15 + j * 49: 25 + j * 49] = +sqrt.(J²[15 + j * 49: 25 + j * 49])
J[26 + j * 49: 44 + j * 49] = -sqrt.(J²[26 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = +sqrt.(J²[45 + j * 49: 49 + j * 49])
#____59
j = 58
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 12 + j * 49] = -sqrt.(J²[06 + j * 49: 12 + j * 49])
J[13 + j * 49: 24 + j * 49] = +sqrt.(J²[13 + j * 49: 24 + j * 49])
J[25 + j * 49: 43 + j * 49] = -sqrt.(J²[25 + j * 49: 43 + j * 49])
J[44 + j * 49: 49 + j * 49] = +sqrt.(J²[44 + j * 49: 49 + j * 49])
#____60
j = 59
J[01 + j * 49: 24 + j * 49] = +sqrt.(J²[01 + j * 49: 24 + j * 49])
J[25 + j * 49: 43 + j * 49] = -sqrt.(J²[25 + j * 49: 43 + j * 49])
J[44 + j * 49: 49 + j * 49] = +sqrt.(J²[44 + j * 49: 49 + j * 49])
end

begin
#____61
j = 60
J[01 + j * 49: 03 + j * 49] = -sqrt.(J²[01 + j * 49: 03 + j * 49])
J[04 + j * 49: 24 + j * 49] = +sqrt.(J²[04 + j * 49: 24 + j * 49])
J[25 + j * 49: 43 + j * 49] = -sqrt.(J²[25 + j * 49: 43 + j * 49])
J[44 + j * 49: 49 + j * 49] = +sqrt.(J²[44 + j * 49: 49 + j * 49])
#____62
j = 61
J[01 + j * 49: 05 + j * 49] = -sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 23 + j * 49] = +sqrt.(J²[06 + j * 49: 23 + j * 49])
J[24 + j * 49: 43 + j * 49] = -sqrt.(J²[24 + j * 49: 43 + j * 49])
J[44 + j * 49: 49 + j * 49] = +sqrt.(J²[44 + j * 49: 49 + j * 49])
#____63
j = 62
J[01 + j * 49: 06 + j * 49] = -sqrt.(J²[01 + j * 49: 06 + j * 49])
J[07 + j * 49: 23 + j * 49] = +sqrt.(J²[07 + j * 49: 23 + j * 49])
J[24 + j * 49: 42 + j * 49] = -sqrt.(J²[24 + j * 49: 42 + j * 49])
J[43 + j * 49: 49 + j * 49] = +sqrt.(J²[43 + j * 49: 49 + j * 49])
#____64
j = 63
J[01 + j * 49: 07 + j * 49] = -sqrt.(J²[01 + j * 49: 07 + j * 49])
J[08 + j * 49: 23 + j * 49] = +sqrt.(J²[08 + j * 49: 23 + j * 49])
J[24 + j * 49: 42 + j * 49] = -sqrt.(J²[24 + j * 49: 42 + j * 49])
J[43 + j * 49: 49 + j * 49] = +sqrt.(J²[43 + j * 49: 49 + j * 49])
#____65
j = 64
J[01 + j * 49: 08 + j * 49] = -sqrt.(J²[01 + j * 49: 08 + j * 49])
J[09 + j * 49: 22 + j * 49] = +sqrt.(J²[09 + j * 49: 22 + j * 49])
J[23 + j * 49: 41 + j * 49] = -sqrt.(J²[23 + j * 49: 41 + j * 49])
J[42 + j * 49: 49 + j * 49] = +sqrt.(J²[42 + j * 49: 49 + j * 49])
#____66
j = 65
J[01 + j * 49: 09 + j * 49] = -sqrt.(J²[01 + j * 49: 09 + j * 49])
J[10 + j * 49: 22 + j * 49] = +sqrt.(J²[10 + j * 49: 22 + j * 49])
J[23 + j * 49: 40 + j * 49] = -sqrt.(J²[23 + j * 49: 40 + j * 49])
J[41 + j * 49: 49 + j * 49] = +sqrt.(J²[41 + j * 49: 49 + j * 49])
#____67
j = 66
J[01 + j * 49: 11 + j * 49] = -sqrt.(J²[01 + j * 49: 11 + j * 49])
J[12 + j * 49: 21 + j * 49] = +sqrt.(J²[12 + j * 49: 21 + j * 49])
J[22 + j * 49: 38 + j * 49] = -sqrt.(J²[22 + j * 49: 38 + j * 49])
J[39 + j * 49: 49 + j * 49] = +sqrt.(J²[39 + j * 49: 49 + j * 49])
#____68
j = 67
J[01 + j * 49: 15 + j * 49] = -sqrt.(J²[01 + j * 49: 15 + j * 49])
J[16 + j * 49: 20 + j * 49] = +sqrt.(J²[16 + j * 49: 20 + j * 49])
J[21 + j * 49: 35 + j * 49] = -sqrt.(J²[21 + j * 49: 35 + j * 49])
J[36 + j * 49: 49 + j * 49] = +sqrt.(J²[36 + j * 49: 49 + j * 49])
#____69
j = 68
J[01 + j * 49: 29 + j * 49] = -sqrt.(J²[01 + j * 49: 29 + j * 49])
J[30 + j * 49: 49 + j * 49] = +sqrt.(J²[30 + j * 49: 49 + j * 49])
#____70
j = 69
J[01 + j * 49: 19 + j * 49] = -sqrt.(J²[01 + j * 49: 19 + j * 49])
J[20 + j * 49: 49 + j * 49] = +sqrt.(J²[20 + j * 49: 49 + j * 49])
end

begin
#____71
j = 70
J[01 + j * 49: 15 + j * 49] = -sqrt.(J²[01 + j * 49: 15 + j * 49])
J[16 + j * 49: 49 + j * 49] = +sqrt.(J²[16 + j * 49: 49 + j * 49])
#____72
j = 71
J[01 + j * 49: 49 + j * 49] = +sqrt.(J²[01 + j * 49: 49 + j * 49])
#____73
j = 72
J[01 + j * 49: 49 + j * 49] = +sqrt.(J²[01 + j * 49: 49 + j * 49])
#____74
j = 73
J[01 + j * 49: 49 + j * 49] = +sqrt.(J²[01 + j * 49: 49 + j * 49])
#____75
j = 74
J[01 + j * 49: 35 + j * 49] = +sqrt.(J²[01 + j * 49: 35 + j * 49])
J[36 + j * 49: 49 + j * 49] = -sqrt.(J²[36 + j * 49: 49 + j * 49])
#____76
j = 75
J[01 + j * 49: 30 + j * 49] = +sqrt.(J²[01 + j * 49: 30 + j * 49])
J[31 + j * 49: 49 + j * 49] = -sqrt.(J²[31 + j * 49: 49 + j * 49])
#____77
j = 76
J[01 + j * 49: 21 + j * 49] = +sqrt.(J²[01 + j * 49: 21 + j * 49])
J[22 + j * 49: 49 + j * 49] = -sqrt.(J²[22 + j * 49: 49 + j * 49])
#____78
j = 77
J[01 + j * 49: 14 + j * 49] = +sqrt.(J²[01 + j * 49: 14 + j * 49])
J[15 + j * 49: 29 + j * 49] = -sqrt.(J²[15 + j * 49: 29 + j * 49])
J[30 + j * 49: 34 + j * 49] = +sqrt.(J²[30 + j * 49: 34 + j * 49])
J[35 + j * 49: 49 + j * 49] = -sqrt.(J²[35 + j * 49: 49 + j * 49])
#____79
j = 78
J[01 + j * 49: 11 + j * 49] = +sqrt.(J²[01 + j * 49: 11 + j * 49])
J[12 + j * 49: 28 + j * 49] = -sqrt.(J²[12 + j * 49: 28 + j * 49])
J[29 + j * 49: 37 + j * 49] = +sqrt.(J²[29 + j * 49: 37 + j * 49])
J[38 + j * 49: 49 + j * 49] = -sqrt.(J²[38 + j * 49: 49 + j * 49])
#____80
j = 79
J[01 + j * 49: 10 + j * 49] = +sqrt.(J²[01 + j * 49: 10 + j * 49])
J[11 + j * 49: 27 + j * 49] = -sqrt.(J²[11 + j * 49: 27 + j * 49])
J[28 + j * 49: 39 + j * 49] = +sqrt.(J²[28 + j * 49: 39 + j * 49])
J[40 + j * 49: 49 + j * 49] = -sqrt.(J²[40 + j * 49: 49 + j * 49])
end

begin
#____81
j = 80
J[01 + j * 49: 09 + j * 49] = +sqrt.(J²[01 + j * 49: 09 + j * 49])
J[10 + j * 49: 26 + j * 49] = -sqrt.(J²[10 + j * 49: 26 + j * 49])
J[27 + j * 49: 41 + j * 49] = +sqrt.(J²[27 + j * 49: 41 + j * 49])
J[42 + j * 49: 49 + j * 49] = -sqrt.(J²[42 + j * 49: 49 + j * 49])
#____82
j = 81
J[01 + j * 49: 08 + j * 49] = +sqrt.(J²[01 + j * 49: 08 + j * 49])
J[09 + j * 49: 26 + j * 49] = -sqrt.(J²[09 + j * 49: 26 + j * 49])
J[27 + j * 49: 42 + j * 49] = +sqrt.(J²[27 + j * 49: 42 + j * 49])
J[43 + j * 49: 49 + j * 49] = -sqrt.(J²[43 + j * 49: 49 + j * 49])
#____83
j = 82
J[01 + j * 49: 07 + j * 49] = +sqrt.(J²[01 + j * 49: 07 + j * 49])
J[08 + j * 49: 26 + j * 49] = -sqrt.(J²[08 + j * 49: 26 + j * 49])
J[27 + j * 49: 43 + j * 49] = +sqrt.(J²[27 + j * 49: 43 + j * 49])
J[44 + j * 49: 49 + j * 49] = -sqrt.(J²[44 + j * 49: 49 + j * 49])
#____84
j = 83
J[01 + j * 49: 06 + j * 49] = +sqrt.(J²[01 + j * 49: 06 + j * 49])
J[07 + j * 49: 25 + j * 49] = -sqrt.(J²[07 + j * 49: 25 + j * 49])
J[26 + j * 49: 44 + j * 49] = +sqrt.(J²[26 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = -sqrt.(J²[45 + j * 49: 49 + j * 49])
#____85
j = 84
J[01 + j * 49: 06 + j * 49] = +sqrt.(J²[01 + j * 49: 06 + j * 49])
J[07 + j * 49: 25 + j * 49] = -sqrt.(J²[07 + j * 49: 25 + j * 49])
J[26 + j * 49: 45 + j * 49] = +sqrt.(J²[26 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____86
j = 85
J[01 + j * 49: 06 + j * 49] = +sqrt.(J²[01 + j * 49: 06 + j * 49])
J[07 + j * 49: 25 + j * 49] = -sqrt.(J²[07 + j * 49: 25 + j * 49])
J[26 + j * 49: 49 + j * 49] = +sqrt.(J²[26 + j * 49: 49 + j * 49])
#____87
j = 86
J[01 + j * 49: 06 + j * 49] = +sqrt.(J²[01 + j * 49: 06 + j * 49])
J[07 + j * 49: 25 + j * 49] = -sqrt.(J²[07 + j * 49: 25 + j * 49])
J[26 + j * 49: 38 + j * 49] = +sqrt.(J²[26 + j * 49: 38 + j * 49])
J[39 + j * 49: 44 + j * 49] = -sqrt.(J²[39 + j * 49: 44 + j * 49])
J[45 + j * 49: 49 + j * 49] = +sqrt.(J²[45 + j * 49: 49 + j * 49])
#____88
j = 87
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 24 + j * 49] = -sqrt.(J²[06 + j * 49: 24 + j * 49])
J[25 + j * 49: 36 + j * 49] = +sqrt.(J²[25 + j * 49: 36 + j * 49])
J[37 + j * 49: 45 + j * 49] = -sqrt.(J²[37 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
#____89
j = 88
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 24 + j * 49] = -sqrt.(J²[06 + j * 49: 24 + j * 49])
J[25 + j * 49: 34 + j * 49] = +sqrt.(J²[25 + j * 49: 34 + j * 49])
J[35 + j * 49: 45 + j * 49] = -sqrt.(J²[35 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
#____90
j = 89
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 24 + j * 49] = -sqrt.(J²[06 + j * 49: 24 + j * 49])
J[25 + j * 49: 33 + j * 49] = +sqrt.(J²[25 + j * 49: 33 + j * 49])
J[34 + j * 49: 45 + j * 49] = -sqrt.(J²[34 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
end

begin
#____91
j = 90
J[01 + j * 49: 05 + j * 49] = +sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 23 + j * 49] = -sqrt.(J²[06 + j * 49: 23 + j * 49])
J[24 + j * 49: 32 + j * 49] = +sqrt.(J²[24 + j * 49: 32 + j * 49])
J[33 + j * 49: 45 + j * 49] = -sqrt.(J²[33 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
#____92
j = 91
J[01 + j * 49: 05 + j * 49] = -sqrt.(J²[01 + j * 49: 05 + j * 49])
J[06 + j * 49: 15 + j * 49] = +sqrt.(J²[06 + j * 49: 15 + j * 49])
J[16 + j * 49: 23 + j * 49] = -sqrt.(J²[16 + j * 49: 23 + j * 49])
J[24 + j * 49: 30 + j * 49] = +sqrt.(J²[24 + j * 49: 30 + j * 49])
J[31 + j * 49: 45 + j * 49] = -sqrt.(J²[31 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = +sqrt.(J²[46 + j * 49: 49 + j * 49])
#____93
j = 92
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 20 + j * 49] = +sqrt.(J²[05 + j * 49: 20 + j * 49])
J[21 + j * 49: 23 + j * 49] = -sqrt.(J²[21 + j * 49: 23 + j * 49])
J[24 + j * 49: 39 + j * 49] = +sqrt.(J²[24 + j * 49: 39 + j * 49])
J[30 + j * 49: 43 + j * 49] = -sqrt.(J²[30 + j * 49: 43 + j * 49])
J[44 + j * 49: 47 + j * 49] = +sqrt.(J²[44 + j * 49: 47 + j * 49])
J[48 + j * 49: 49 + j * 49] = -sqrt.(J²[48 + j * 49: 49 + j * 49])
#____94
j = 93
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 22 + j * 49] = +sqrt.(J²[05 + j * 49: 22 + j * 49])
J[23 + j * 49: 38 + j * 49] = -sqrt.(J²[23 + j * 49: 38 + j * 49])
J[39 + j * 49: 46 + j * 49] = +sqrt.(J²[39 + j * 49: 46 + j * 49])
J[47 + j * 49: 49 + j * 49] = -sqrt.(J²[47 + j * 49: 49 + j * 49])
#____95
j = 94
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 22 + j * 49] = +sqrt.(J²[05 + j * 49: 22 + j * 49])
J[23 + j * 49: 31 + j * 49] = -sqrt.(J²[23 + j * 49: 31 + j * 49])
J[32 + j * 49: 45 + j * 49] = +sqrt.(J²[32 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____96
j = 95
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 21 + j * 49] = +sqrt.(J²[05 + j * 49: 21 + j * 49])
J[22 + j * 49: 30 + j * 49] = -sqrt.(J²[22 + j * 49: 30 + j * 49])
J[31 + j * 49: 45 + j * 49] = +sqrt.(J²[31 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
#____97
j = 96
J[01 + j * 49: 04 + j * 49] = -sqrt.(J²[01 + j * 49: 04 + j * 49])
J[05 + j * 49: 20 + j * 49] = +sqrt.(J²[05 + j * 49: 20 + j * 49])
J[21 + j * 49: 29 + j * 49] = -sqrt.(J²[21 + j * 49: 29 + j * 49])
J[30 + j * 49: 45 + j * 49] = +sqrt.(J²[30 + j * 49: 45 + j * 49])
J[46 + j * 49: 49 + j * 49] = -sqrt.(J²[46 + j * 49: 49 + j * 49])
end

begin
	J²_mat = zeros((49,Int(size(J²,1)/49)))
	J_mat = zeros((49,Int(size(J²,1)/49)))
	logJ_mat = zeros((49,Int(size(J²,1)/49)))
	θ = zeros((Int(size(J²,1)/49)))
	ϕ = zeros((49))
	for bb in 1:Int(size(J²,1)/49)
		J²_mat[:, bb] = J²[49*(bb-1)+1: 49*bb]
		J_mat[:, bb] = J[49*(bb-1)+1: 49*bb]
		logJ_mat[:, bb] = logJ²[49*(bb-1)+1: 49*bb]
		θ[bb] = para[49*bb,1]
	end
	for aa in 1:49
		ϕ[aa] = para[aa,2]
	end
end

anim2 = @animate for ii in 1:97
    plot!(scatter(ϕ, J_mat[:,ii], ylims = (-0.5,0.5)),
	xlabel= L"$\phi$", ylabel=L"J",
	xticks=([0, π, π/2], [L"0", L"\pi", L"\pi/2"]),
	title=L"J " * L"$(\theta = $" *string(round(θ[ii], digits=2)) * L")", legend = false)
    annotations = [(ϕ[i], J[i+(ii-1)*49], text(string(i), :left, 12)) for i in 1:length(ϕ)]
    annotate!(annotations)
end
gif(anim2, "anim2.gif", fps = 20)


anim3 = @animate for ii in 1:49
    plot!(scatter(θ, J_mat[ii, :], ylims = (-0.5,0.5)),
	xlabel= L"$\theta$", ylabel=L"J",
	xticks=([0, -π, π, -π/2, π/2], ["0", L"-\pi", L"\pi", L"-\pi/2",L"\pi/2"]),
	title=L"J " * L"$(\phi = $" *string(round(ϕ[ii], digits=2)) *L")", legend = false)
    annotations = [(θ[i], J_mat[ii, i], text(string(i), :left, 12)) for i in 1:length(θ)]
    annotate!(annotations)
end

gif(anim3, "anim3.gif", fps = 20)


anim4 = @animate for ii in 1:97
    plot!(scatter(ϕ, J²_mat[:,ii], ylims = (-0.01,0.25)),
	xlabel= L"$\phi$", ylabel=L"J^2",
	xticks=([0, π, π/2], [L"0", L"\pi", L"\pi/2"]),
	title=L"J^2 " * L"$(\theta = $" *string(round(θ[ii], digits=2)) * L")", legend = false)
    annotations = [(ϕ[i], J²[i+(ii-1)*49], text(string(i), :left, 12)) for i in 1:length(ϕ)]
    annotate!(annotations)
end
gif(anim4, "anim4.gif", fps = 30)

anim5 = @animate for ii in 1:49
    plot!(scatter(θ, J²_mat[ii, :], ylims = (-0.01,0.25)),
	xlabel= L"$\theta$", ylabel=L"J^2",
	xticks=([0, -π, π, -π/2, π/2], ["0", L"-\pi", L"\pi", L"-\pi/2",L"\pi/2"]),
	title=L"J^2 " * L"$(\phi = $" *string(round(ϕ[ii], digits=2)) *L")", legend = false)
    annotations = [(θ[i], J²_mat[ii, i], text(string(i), :left, 12)) for i in 1:length(θ)]
    annotate!(annotations)
end
gif(anim5, "anim5.gif", fps = 30)

begin
plot(θ, ϕ,
	J_mat,
	st=:surface,
	xlabel= L"$\phi$", ylabel=L"$\theta$", zlabel=L"J",
	yticks=([0, π, π/2], [L"0", L"\pi", L"\pi/2"]),
	xticks=([0, -π, π, -π/2, π/2], ["0", L"-\pi", L"\pi", L"-\pi/2",L"\pi/2"]),
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"])
	)
savefig("J_surface.png")
end

begin
heatmap(ϕ[1:49],
	θ[1:97], J_mat[:,:],
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
	xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	xticks=([0, π, π/2], [L"0", L"\pi", L"\pi/2"]),
	yticks=([0, -π, π, -π/2, π/2], ["0", L"-\pi", L"\pi", L"-\pi/2",L"\pi/2"]),
	#aspectratio=:equal,
	size=(700, 700),
	title="Transfer Integrals" )
	savefig("J_heatmap.png")
end

begin
plot(θ,ϕ,
	J²_mat,
	st=:surface,
	xlabel= L"$\phi$", ylabel=L"$\theta$", zlabel=L"J",
	yticks=([0, π, π/2], [L"0", L"\pi", L"\pi/2"]),
	xticks=([0, -π, π, -π/2, π/2], ["0", L"-\pi", L"\pi", L"-\pi/2",L"\pi/2"]),
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"])
	)
	savefig("J^2_surface.png")
end

begin
heatmap(ϕ[1:49],
	θ[1:97], J²_mat[:,:],
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
	xlabel= L"$\phi$", ylabel=L"$\theta$",
	xticks=([0, π, π/2], [L"0", L"\pi", L"\pi/2"]),
	yticks=([0, -π, π, -π/2, π/2], ["0", L"-\pi", L"\pi", L"-\pi/2",L"\pi/2"]),
	#aspectratio=:equal,
	size=(700, 700),
	title="Squared Transfer Integrals" )
	savefig("J^2_heatmap.png")
end

begin
	function rot_x(θ)
		return [[1.0, 0.0, 0.0] [0.0, cos(θ), sin(θ)] [0.0, -sin(θ), cos(θ)]]
	end
	function rot_y(θ)
		return [[cos(θ), 0.0, -sin(θ)] [0.0, 1.0, 0.0] [sin(θ), 0.0, cos(θ)]]
	end
	function rot_z(θ)
		return [[cos(θ), sin(θ), 0.0] [-sin(θ), cos(θ), 0.0] [0.0, -sin(θ), 1.0]]
	end
end;

function construct_dimer(namefile, savefile, para, rot_ang_y, J)
		mol1 = read_extxyz(namefile)[1] 	#Read Molecular Structure from file.extxyz
		ndata = size(para, 1) 				#Collect number of data from TranInt out.dat
		natoms = size(mol1.M)[1] 			#Collect number of atoms in molecules
		Set_Mol2 = Vector{Atoms{Float64}}([])						#Set of 2nd Molecules
		for jj in 1:ndata
			mol2 = read_extxyz(namefile)[1]
			#println("Dimer: "*string(jj))
			for ii in 1:natoms
				mol2.X[ii] = rot_y(rot_ang_y) * mol1.X[ii] + para[jj, 4:6] 		#Rotation + Translation
				#mol2.X[ii] = rot_z(para[jj, 2])*rot_y(para[jj, 1]) * mol1.X[ii] + para[jj, 4:6] 		#Rotation + Translation
			#println(mol2.X[ii])
				set_data!(mol2, "energy", J[jj])
			end
			push!(Set_Mol2, mol2) 	
			#println(Set_Mol2)
			#println("_________________________________________________________")
		end
		#return Set_Mol2
		write_extxyz(savefile, Set_Mol2)
end

construct_dimer("naphthalene.extxyz",
"dimer_naphthalene.extxyz", para, 2*pi/3, J);

begin
	function find_element(Molecule, Z)
		return findall(x->x==Z, Molecule.Z)
	end
	mol1 = read_extxyz("naphthalene.extxyz")[1]
	num_C = find_element(mol1, 6)
	shared_C = [3, 4];
end;

function Centroid_bead(Molecule, shared_C, num_C)
	cent_Vec = [] 					#collect all bead centroids
	
	ring_C = zeros(Int64,(2,6))		#collect Carbon atoms (Shared + Unshared)
	ring_C[1, 1:6] = num_C[1:6]
	ring_C[2, 1:2] = shared_C
	ring_C[2, 3:6] = num_C[7:10]
	
	for jj in 1:2
		vec = zeros((3))
		bead_mass = 0
		for ii in 1:6
			vec += Molecule.X[ring_C[jj, ii]]*Molecule.M[ring_C[jj, ii]]
			bead_mass += Molecule.M[ring_C[jj, ii]]
		end
		#println(vec)
		#println(bead_mass)
		#println(vec/bead_mass)
		#println(vec/78.11)
		push!(cent_Vec, vec/bead_mass) #collect centroid positions with the CM of the beads
	end
	return cent_Vec
end


function construct_all_bead(namefile1, namefile2, savefile, shared_C, num_C, atomic_Z)
	
	mol1 = read_extxyz(namefile1)[1] 							#first mole
	Set_Mol2 =  read_extxyz(namefile2)							#set of second mole
	
	x1 = Centroid_bead(mol1, shared_C, num_C)				#Centroid first
	bead_mass = fill(sum(mol1.M)/2, 4)						#Set bead mass
	
	Large_data = Vector{Atoms{Float64}}([])					#write in extxyz

	#types = fill(0, 4)
	types = fill(atomic_Z, 4)
	
	for nn in 1:length(Set_Mol2)
		mol2 = Set_Mol2[nn]
		x2_n = Centroid_bead(mol2, shared_C, num_C)

		Dimer = Atoms(X = vcat(x2_n, x1),
			M = bead_mass,
			Z = types,
			cell = mol1.cell*100,
			pbc = fill(false, 3)			
		)

		set_data!(Dimer, "energy", 
			get_data(mol2, "energy")
		)
		
		push!(Large_data, Dimer)
	end
	write_extxyz(savefile, Large_data)
end

construct_all_bead("naphthalene.extxyz", 
	"dimer_naphthalene.extxyz", 
	"Large_bead.extxyz", 
	shared_C, 
	num_C,
	14
)


function construct_all_dimer(namefile1, namefile2, savefile)
	mol1 = read_extxyz(namefile1)[1] 			#Read Molecular Structure from 													file.extxyz
	Set_Mol2 =  read_extxyz(namefile2)	
	Set_dimer = Vector{Atoms{Float64}}([])		#Set of dimer Molecules
	
	for nn in 1:length(Set_Mol2)
		mol2 = Set_Mol2[nn]
		#println("Dimer: "*string(jj))
		dimer = Atoms(X = vcat(mol2.X, mol1.X),
			M = vcat(mol2.M, mol1.M),
			Z = vcat(mol2.Z, mol1.Z),
			cell = mol1.cell*100,
			pbc = fill(false, 3)
		)
		
		set_data!(dimer, "energy", 
			get_data(mol2, "energy")
		)
		
		push!(Set_dimer, dimer) 	
		#println(Set_Mol2)
		#println("_________________________________________________________")
	end
	#return Set_Mol2
	write_extxyz(savefile, Set_dimer)
end

construct_all_dimer("naphthalene.extxyz", 
	"dimer_naphthalene.extxyz", 
	"all_CH_dimers.extxyz"
)

function train_test(all_structure_file, num_train, train_file, test_file)
	s_dim =  read_extxyz(all_structure_file)
	
	num_test = length(s_dim) - num_train

	trian_in = ones((1))

	for nn in 1:num_train
		ind = rand(1:length(s_dim))
		if !(ind in trian_in)
			append!(trian_in, ind)
		end
	end
	println(length(trian_in))
	train_index = Int.(trian_in)
	test_index = setdiff(collect(1:length(s_dim)), train_index)
	
	train_dat = Vector{Atoms{Float64}}([])	
	test_dat = Vector{Atoms{Float64}}([])

	for ii in train_index
		confi = s_dim[ii]
		push!(train_dat, confi)
	end
	write_extxyz(train_file, train_dat)
	for jj in test_index
		confi = s_dim[jj]
		push!(test_dat, confi)
	end
	write_extxyz(test_file, test_dat)
	
end

train_test("all_CH_dimers.extxyz", 2000, "train_all.extxyz", "test_all.extxyz")
train_test("Large_bead.extxyz", 2000, "train_bead.extxyz", "test_bead.extxyz")

function construct_Carbon(oldfile, newfile, J)
	Set_Mol2 =  read_extxyz(oldfile)	
	Set_dimer = Vector{Atoms{Float64}}([])		#Set of dimer Molecules
	
	for nn in 1:length(Set_Mol2)
		mol2 = Set_Mol2[nn]
		
		set_data!(mol2, "energy", 
			J[nn]
		)
		push!(Set_dimer, mol2) 	
		#println(Set_Mol2)
		#println("_________________________________________________________")
	end
	#return Set_Mol2
	write_extxyz(newfile, Set_dimer)
end

begin
	Dimer_full = read_extxyz("dimer_naphthalene.extxyz")
	    
	mol = Dimer_full[1]
	index = findall(x -> x == 6, mol.Z)
	num = length(index)
	    
	if num == 0
	    println("Warning: No atoms found with atomic number $Z.")
	    return
	end
	    
	Mass = fill(mol.M[1], num)
	Z_atomic = fill(mol.Z[1], num)
	cell = mol.cell
	    
	Large_data = Vector{Atoms{Float64}}([])
end

begin
	for jj in 1:length(Dimer_full)
		mol1 = Dimer_full[jj]
        X1 = [] 
		
        for ii in index
            append!(X1, mol1.X[ii])
        end
	dimer = Atoms(X = X1,
	                M = Mass,
	                Z = Z_atomic,
	                cell = cell,
	                pbc = fill(false, 3))

	set_data!(dimer, "energy", 
			get_data(mol1, "energy")
		)
		
	push!(Large_data, dimer)
	end
end
write_extxyz("Large_dimers_carbon.extxyz", Large_data)

function construct_all_carbon(namefile1, namefile2, savefile)
	mol1 = read_extxyz(namefile1)[1]
	Set_Mol2 =  read_extxyz(namefile2)	
	Set_dimer = Vector{Atoms{Float64}}([])

	X1 = [] 
    for ii in index
        append!(X1, mol1.X[ii])
    end

	mol = Atoms(X = X1,
				M = Mass,
				Z = Z_atomic,
				cell = cell,
				pbc = fill(false, 3))

	for nn in 1:length(Set_Mol2)
		mol2 = Set_Mol2[nn]
		#println("Dimer: "*string(jj))
		dimer = Atoms(X = vcat(mol2.X, mol.X),
			M = vcat(mol2.M, mol.M),
			Z = vcat(mol2.Z, mol.Z),
			cell = mol.cell*100,
			pbc = fill(false, 3)
		)
		
		set_data!(dimer, "energy", 
			get_data(mol2, "energy")
		)
		
		push!(Set_dimer, dimer) 	
		#println(Set_Mol2)
		#println("_________________________________________________________")
	end
	#return Set_Mol2
	write_extxyz(savefile, Set_dimer)
end
construct_all_carbon("naphthalene.extxyz", 
		"Large_dimers_carbon.extxyz",
		"dimers_carbon.extxyz"
		)

train_test("dimers_carbon.extxyz", 2000, "train_carbon.extxyz", "test_carbon.extxyz")

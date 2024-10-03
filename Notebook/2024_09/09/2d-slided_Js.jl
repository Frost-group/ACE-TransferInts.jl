
using CSV
using DataFrames
using DelimitedFiles

function read_J(namefile)
	data_J = (CSV.File(namefile, delim=" ")|> Tables.matrix)
	para, J, logJ = data_J[:, 1:6], data_J[:, 7], data_J[:, 8]
	return para, J, logJ
end

#Reading parameters the separation of 2 molecular centers (θ, ϕ, r, x, y, z), transfer integrals Js and log(Js),
#where (θ, ϕ, r) is the spherical coordinates transformed from (x, y, z).
para, J, logJ = read_J("out.dat")

#Creating Vectors and Matrices to contain parameters, Js and log(Js)
J_mat = zeros((49,Int(size(J,1)/49)))
logJ_mat = zeros((49,Int(size(J,1)/49)))
θ = zeros((Int(size(J,1)/49)))
ϕ = zeros((49))


for bb in 1:Int(size(J,1)/49)
	J_mat[:, bb] = J[49*(bb-1)+1: 49*bb]	
	logJ_mat[:, bb] = logJ[49*(bb-1)+1: 49*bb]
	θ[bb] = para[49*bb,1] #running from -π to π
end

for aa in 1:49
	ϕ[aa] = para[aa,2] #running from 0 to π
end

using Plots
using LaTeXStrings

heatmap(ϕ[1:49],
	θ[1:97], J_mat[:,:],
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
	xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	#aspectratio=:equal,
	size=(700, 700),
	title="Transfer Integrals" )


plot(θ,ϕ,
    logJ_mat,
    st=:surface,
    xlabel= L"$\phi$",
    ylabel=L"$\theta$",
    zlabel=L"J",
    c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]))


θ_i = [10,30,50,70]
ϕ_i = [10, 30, 40]
Color = ["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]

begin
	plt1 = plot()
	for n = 1:4
		θ_j = round(θ[θ_i[n]], digits = 2)	
		plot!(ϕ, logJ_mat[:,θ_i[n]], lw = 2, c = Color[n]
			, label = L"\theta ="*" $θ_j"
			, xlabel = L"\phi", ylabel = L"log(J)"
		)
	end
	plt1
	plt2 = plot()
	for m = 1:3
		ϕ_j = round(ϕ[ϕ_i[m]], digits = 2)	
		plot!(θ, logJ_mat[ϕ_i[m],:], lw = 2, c = Color[2*m-1]
				, label = L"\theta ="*" $ϕ_j"
				, xlabel = L"\theta", ylabel = L"log(J)"
		)
	end
	plt2
	plot!(plt1, plt2, layout=(2,1), size=(900,700))
	savefig("1d-sliced_logJ.png")
end


plot(θ,ϕ,
    logJ_mat,
    st=:surface,
    xlabel= L"$\phi$",
    ylabel=L"$\theta$",
    zlabel=L"J",
    c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]))

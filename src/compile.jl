using MATH50003NumericalAnalysis


####
# Notes
####

compilenotes("I.1.RectangularRule")
compilenotes("I.2.DividedDifferences")
compilenotes("I.3.DualNumbers")
compilenotes("I.4.NewtonMethod")

compilenotes("II.1.Reals")
compilenotes("II.2.Arithmetic")
compilenotes("II.3.Intervals")

compilenotes("III.1.StructuredMatrices")
# compilenotes("III.2.DifferentialEquations")
# compilenotes("III.3.Cholesky")
# compilenotes("III.4.Regression")
# compilenotes("III.5.OrthogonalMatrices")
# compilenotes("III.6.QR")

# compilenotes("IV.1.Fourier")
# compilenotes("IV.2.DFT")
# compilenotes("IV.3.OrthogonalPolynomials")
# compilenotes("IV.4.ClassicalOPs")
# compilenotes("IV.5.GaussianQuadrature")


compilenotes("A.Asymptotics")
compilenotes("A.Integers")
compilenotes("A.Permutations")


####
# Sheets
####

for k = 1:3
    compilesheet(k)
end


for k = 1:2
    compilesheetsolution(k)
end

####
# Labs
####
for k = 1:3
    compilelab(k)
    compilelabdemo(k)
end




compilelabsolution(1)
compilelabsolution(2)
# compilelabsolution(3)
# compilelabsolution(4)
# compilelabsolution(5)
# compilelabsolution(6)
# compilelabsolution(7)
# compilelabsolution(8)



#####
# extras
#####

using Weave

nkwds = (out_path="notes/", jupyter_path="$(homedir())/.julia/conda/3/x86_64/bin/jupyter", nbconvert_options="--allow-errors")
notebook("src/notes/A.Julia.jmd"; nkwds...)


###
# exams
###

import Literate

function compileexam(str)
    write("exams/$str.jl", replace(replace(read("src/exams/$(str)s.jl", String), r"## SOLUTION(.*?)## END"s => "")))
    Literate.notebook("exams/$str.jl", "exams/"; execute=false)
end

function compileexamsolution(str)
    Literate.notebook("src/exams/$(str)s.jl", "exams/")
end

compileexam("mockexam")
compileexamsolution("mockexam")
function sol = evaluate_helmholtz(data, mesh, sample)

if mesh.homogeneous_flag
    sol = evaluate_helmholtz_constantk(data, mesh, sample);
else
    sol = evaluate_helmholtz_variablek(data, mesh, sample);
end








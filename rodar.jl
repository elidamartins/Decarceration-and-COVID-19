# rodar_1000_vezes.jl

# Definir o número de vezes que o script será executado
n_vezes = 100

# Loop para executar o script principal várias vezes
for i in 1:n_vezes
    # Comando para chamar o script principal usando `include`
    include("simulacao_presidio48.jl")
end
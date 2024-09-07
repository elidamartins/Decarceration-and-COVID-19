using Random
using CSV
using DataFrames
using Dates
using DelimitedFiles

date_format = "yyyy-mm-dd_HH-MM-SS"
time = Dates.format(now(), date_format)

code2state_name = Dict(0 => "suscetivel",
                       1 => "infectado",
                       2 => "exposto",
                       3 => "vacinado",
                       4 => "grave",
                       5 => "critico",
                       6 => "morto_covid",
                       7 => "saida",
                       8 => "recuperado");

state_name2code = Dict(value => key for (key, value) in code2state_name)

mutable struct Preso
    cela
    idade
    estado
    dia_alteracao_estado
    proximo_estado
    alterar_apos
    grupo
    assintomatico # bool
    dias_vividos # Adicionando campo dias_vividos
end

mutable struct Prisao
    presos
    tempo_de_sol # em horas
end

mutable struct Mortos
    presos
end

mutable struct Hospital
    presos
    tempo_desde_chegada
end

mutable struct Simulacao
    tempo_de_simulacao::Int
    n_presos::Int
    prisao::Prisao
    hospital::Hospital
    mortos::Mortos
    tx_entrada::Int
    tx_saida::Int
    tx_entrada_inicial::Int
    tx_saida_inicial::Int
    mortalidade_hospitalar::Float64
    tx_infectado_grave::Vector{Float64}
    tx_grave_critico::Float64
    vacinacao_por_dia::Int
    dia::Int
    probabilidade_contagio::Float64
    n_celas::Int
    dias_de_vacinacao::Vector{Int}
    distribuicao_etaria::Vector{Float64}
    n_grupos::Int
    tx_assintomaticos::Float64
    taxa_encontros_celas::Float64
    p_vacinacao::Float64
    dias_recuperados_para_sucetiveis::Vector{Float64}
    dias_vacinados_para_sucetiveis::Vector{Float64}
    dias_expostos_para_infectados::Vector{Float64}
    dias_infectados_para_graves::Vector{Float64}
    dias_graves_para_criticos::Vector{Float64}
    dias_criticos_para_mortos::Vector{Float64}
    dias_infectados_para_recuperados::Vector{Float64}
    dias_graves_para_recuperados::Vector{Float64}
    dias_criticos_para_recuperados::Vector{Float64}
    custo_total::Float64
    total_pessoas::Int
end

# Definir os custos diários para cada estado
estado_custos = Dict(
    state_name2code["suscetivel"] => 90.46,
    state_name2code["infectado"] => 107.61,
    state_name2code["exposto"] => 90.46,
    state_name2code["vacinado"] => 90.46,
    state_name2code["grave"] => 537.63,
    state_name2code["critico"] => 2613.04,
    state_name2code["morto_covid"] => 0,
    state_name2code["saida"] => 0,
    state_name2code["recuperado"] => 90.46
)

# Definir o custo adicional para cada entrada no hospital
custo_entrada_hospital = 273.65  # ajuste conforme necessário
custo_entrada_vacinacao = 72.13  # ajuste conforme necessário

expectativa_vida_dias = 1092


function calcular_custo_diario(simulacao, identificador)
    custo_total = 0
    for preso in simulacao.prisao.presos
        custo_total += estado_custos[preso.estado]
    end
    for preso in simulacao.hospital.presos
        custo_total += estado_custos[preso.estado]
    end
    return custo_total
end
#=
function registrar_custos_diarios(custo, identificador)
    open("custos_diarios_$identificador.txt", "a") do file
        write(file, string(custo) * "\n")
    end
end
=#
#=
function registrar_custo_total(custo_total, identificador)
    open("custo_total_$identificador.txt", "w") do file
        write(file, string(custo_total) * "\n")
    end
end
=#
function aplicar_taxa_desconto(custo_total, semanas, taxa_anual)
    anos = semanas / 52
    custo_descontado = custo_total / (1 + taxa_anual)^anos
    return custo_descontado
end
#=
function registrar_custo_total_descontado(custo_descontado, identificador)
    open("custo_total_descontado_$identificador.txt", "w") do file
        write(file, string(custo_descontado) * "\n")
    end
end
=#
function calcular_custo_descontado_por_individuo(custo_descontado, simulacao)
    custo_descontado_por_individuo = custo_descontado / simulacao.total_pessoas
    return custo_descontado_por_individuo
end

#=
function registrar_custo_descontado_por_individuo(custo_descontado_por_individuo, identificador)
    open("custo_descontado_por_individuo_$identificador.txt", "w") do file
        write(file, string(custo_descontado_por_individuo) * "\n")
    end
end
=#

function sorteio(pmf, espaco_amostral) 
    acumulada = accumulate(+,pmf)
    moeda = rand()
    resultado  = nothing
    i = 1
    while resultado == nothing 
        if (moeda <= acumulada[i])
            resultado = espaco_amostral[i]
        else 
            i = i+1
        end
    end
    return resultado
end

function instanciar_Preso(cela, distribuicao_etaria, n_grupos)
    espaco_amostral = [1, 2, 3, 4, 5, 6, 7]
    idades_por_faixa = [18 24; 25 29; 30 35; 36 45; 46 60; 61 70; 71 80]
    faixa = sorteio(distribuicao_etaria, espaco_amostral)
    idade = rand(idades_por_faixa[faixa, 1]:idades_por_faixa[faixa, 2])
    grupo = Int(rand(1:n_grupos))
    return Preso(cela, idade, state_name2code["suscetivel"], 1, nothing, 0, grupo, false, 0)
end

function instanciar_prisao(tempo_sol, n_celas, total_de_presos, distribuicao_etaria, n_grupos)
    presos = []
    
    presos_por_cela = total_de_presos ÷ n_celas
    resto_de_presos = total_de_presos % n_celas

    for i in 1:n_celas
        if i == n_celas
            for j in 1:(presos_por_cela + resto_de_presos)
                preso = instanciar_Preso(i, distribuicao_etaria, n_grupos)
                push!(presos, preso)
            end
        else
            for j in 1:presos_por_cela
                preso = instanciar_Preso(i, distribuicao_etaria, n_grupos)
                push!(presos, preso)
            end
        end
    end

    return Prisao(presos, tempo_sol)
end

function gerar_encontros(N, taxa_encontros)
    N_nao_nulos = Int(round(taxa_encontros * (N^2 - N)))

    if N_nao_nulos == 0
        encontros = Int64[]
    else
        encontros = rand(1:N, N_nao_nulos, 2) # preso i encontrou preso j
    end
    return encontros
end

function simular_banho_de_sol(simulacao, identificador)
    taxa_encontros = simulacao.prisao.tempo_de_sol / 15

    for idx_grupo in 1:simulacao.n_grupos
        grupo = cria_grupo(simulacao, idx_grupo)

        # Intra-grupo
        encontros = gerar_encontros(length(grupo), taxa_encontros)
        for linha in 1:size(encontros, 1)
            i = encontros[linha, 1]
            j = encontros[linha, 2]

            contagio, posicao = houve_contagio(grupo, i, j, simulacao.probabilidade_contagio)
            if contagio
                grupo[posicao].estado = state_name2code["exposto"]
                grupo[posicao].proximo_estado = state_name2code["infectado"] 
                espaco_amostral = collect(simulacao.dias_expostos_para_infectados[1]:1:simulacao.dias_expostos_para_infectados[3])
                sdt1emeio = length(espaco_amostral)/3
                pmf = exp.(-((espaco_amostral.-simulacao.dias_expostos_para_infectados[2]).^2)/sdt1emeio)
                pmf = pmf/sum(pmf)
                grupo[posicao].alterar_apos = sorteio(pmf, espaco_amostral) #vai para o estado infectado depois de 4 ou 5 dias
                grupo[posicao].dia_alteracao_estado = simulacao.dia
            end
        end
    end

    # Inter-grupo
    taxa_encontros = simulacao.prisao.tempo_de_sol / 1999
    encontros = gerar_encontros(length(simulacao.prisao.presos), taxa_encontros)
    for linha in 1:size(encontros, 1)
        i = encontros[linha, 1]
        j = encontros[linha, 2]

        contagio, posicao = houve_contagio(simulacao.prisao.presos, i, j, simulacao.probabilidade_contagio)

        if contagio
            simulacao.prisao.presos[posicao].estado = state_name2code["exposto"]
            simulacao.prisao.presos[posicao].proximo_estado = state_name2code["infectado"] 
            espaco_amostral = collect(simulacao.dias_expostos_para_infectados[1]:1:simulacao.dias_expostos_para_infectados[3])
            sdt1emeio = length(espaco_amostral)/3
            pmf = exp.(-((espaco_amostral.-simulacao.dias_expostos_para_infectados[2]).^2)/sdt1emeio)
            pmf = pmf/sum(pmf)
            simulacao.prisao.presos[posicao].alterar_apos = sorteio(pmf, espaco_amostral) #vai para o estado infectado depois de 4 ou 5 dias
            simulacao.prisao.presos[posicao].dia_alteracao_estado = simulacao.dia
        end
    end
end

function cria_cela(simulacao, idx_cela)
    cela = []

    for i in 1:length(simulacao.prisao.presos)
        if simulacao.prisao.presos[i].cela == idx_cela
            push!(cela, simulacao.prisao.presos[i])
        end
    end
    
    return cela
end

function cria_grupo(simulacao, idx_grupo)
    grupo = []

    for i in 1:length(simulacao.prisao.presos)
        if simulacao.prisao.presos[i].grupo == idx_grupo
            push!(grupo, simulacao.prisao.presos[i])
        end
    end
    
    return grupo
end

function houve_contagio(presos, preso_i, preso_j, probabilidade)
    contagio = (
        (presos[preso_i].estado == 1 && presos[preso_j].estado == 0) ||
        (presos[preso_i].estado == 0 && presos[preso_j].estado == 1) 
    ) && rand() < probabilidade;

    if contagio
        if presos[preso_i].estado == 1
            return contagio, preso_j
        else
            return contagio, preso_i
        end
    else
        return contagio, nothing
    end
end

function simular_celas(simulacao, identificador)
    taxa_de_encontros = simulacao.taxa_encontros_celas

    for idx_cela in 1:simulacao.n_celas
        cela = cria_cela(simulacao, idx_cela)
        encontros = gerar_encontros(length(cela), taxa_de_encontros)
        for linha in 1:size(encontros,1)
            i = encontros[linha, 1]
            j = encontros[linha, 2]

            contagio, posicao = houve_contagio(cela, i, j, simulacao.probabilidade_contagio)
            if contagio
                cela[posicao].estado = state_name2code["exposto"]
                cela[posicao].proximo_estado = state_name2code["infectado"]
                espaco_amostral = collect(simulacao.dias_expostos_para_infectados[1]:1:simulacao.dias_expostos_para_infectados[3])
                sdt1emeio = length(espaco_amostral)/3
                pmf = exp.(-((espaco_amostral.-simulacao.dias_expostos_para_infectados[2]).^2)/sdt1emeio)
                pmf = pmf/sum(pmf)
                cela[posicao].alterar_apos = sorteio(pmf, espaco_amostral) #vai para o estado infectado depois de 4 ou 5 dias
                cela[posicao].dia_alteracao_estado = simulacao.dia
            end
        end
    end
end

#=
function registrar_suscetiveis(simulacao, identificador)
    # Conta o número de indivíduos suscetíveis
    num_suscetiveis = count(p -> p.estado == state_name2code["suscetivel"], simulacao.prisao.presos)
    
    # Registra o número de suscetíveis no arquivo
    open("suscetiveis_por_dia_$identificador.txt", "a") do file
        write(file, string(num_suscetiveis) * "\n")
    end
end
=#
#=
function registrar_expostos(simulacao, identificador)
    # Conta o número de indivíduos expostos
    num_expostos = count(p -> p.estado == state_name2code["exposto"], simulacao.prisao.presos)
    
    # Registra o número de expostos no arquivo
    open("expostos_por_dia_$identificador.txt", "a") do file
        write(file, string(num_expostos) * "\n")
    end
end
=#
#=
function registrar_infectados(simulacao, identificador)
    # Conta o número de indivíduos infectados
    num_infectados = count(p -> p.estado == state_name2code["infectado"], simulacao.prisao.presos)
    
    # Registra o número de infectados no arquivo
    open("infectados_por_dia_$identificador.txt", "a") do file
        write(file, string(num_infectados) * "\n")
    end
end
=#
function expostos_para_infectados(simulacao, identificador)
    idx_infectar = findall(a->a.proximo_estado==state_name2code["infectado"] && (a.dia_alteracao_estado+a.alterar_apos == simulacao.dia), simulacao.prisao.presos) 
    open("casos_por_dia_$identificador.txt", "a") do file
        write(file, string(length(idx_infectar))*"\n")
    end
    faixas = [18, 25,46,61,70]
    for i in idx_infectar
        idade = simulacao.prisao.presos[i].idade

        #escolhendo a probabilidade de agravar dependendo da idade
        j = 1
        faixa = nothing
        while faixa == nothing
            if idade < faixas[j+1]
                faixa = j 
            else 
                j = j+1
            end 
            if j == length(faixas)
                faixa = j
            end
        end   
        simulacao.prisao.presos[i].estado = state_name2code["infectado"];
        simulacao.prisao.presos[i].dia_alteracao_estado = simulacao.dia;
        #definindo o próximo estado
        p_hosp = simulacao.tx_infectado_grave[faixa]*(1-simulacao.tx_assintomaticos)
        espaco_amostral = ["recuperado", "grave"]
        pmf = [1-p_hosp, p_hosp]
        destino = sorteio(pmf, espaco_amostral)
        simulacao.prisao.presos[i].proximo_estado = state_name2code[destino]

        #definindo apos quantos dias haverá mudança de estado
        if (destino == "recuperado")
            espaco_amostral = collect(simulacao.dias_infectados_para_recuperados[1]:1:simulacao.dias_infectados_para_recuperados[3])
            sdt1emeio = length(espaco_amostral)/3
            pmf = exp.(-((espaco_amostral.-simulacao.dias_infectados_para_recuperados[2]).^2)/sdt1emeio)
            pmf = pmf/sum(pmf)
            dias_para_mudar = sorteio(pmf, espaco_amostral)
            simulacao.prisao.presos[i].alterar_apos = dias_para_mudar
        else
            espaco_amostral = collect(simulacao.dias_infectados_para_graves[1]:1:simulacao.dias_infectados_para_graves[3])
            sdt1emeio = length(espaco_amostral)/3
            pmf = exp.(-((espaco_amostral.-simulacao.dias_infectados_para_graves[2]).^2)/sdt1emeio)
            pmf = pmf/sum(pmf)
            simulacao.prisao.presos[i].alterar_apos = sorteio(pmf, espaco_amostral)
        end
    end
end

function infectados_para_graves(simulacao, identificador)
    idx_agravar = findall(a -> (a.proximo_estado == state_name2code["grave"]) && (a.dia_alteracao_estado + a.alterar_apos == simulacao.dia), simulacao.prisao.presos)
    #open("entradas_hospital_por_dia_$identificador.txt", "a") do file
        #write(file, string(length(idx_agravar)) * "\n")
    #end
    for i in idx_agravar
        preso = simulacao.prisao.presos[i]
        preso.dia_alteracao_estado = simulacao.dia
        preso.estado = state_name2code["grave"]
        espaco_amostral = ["recuperado", "critico", "morto_covid"]
        pmf = [1 - simulacao.tx_grave_critico - simulacao.mortalidade_hospitalar, simulacao.tx_grave_critico, simulacao.mortalidade_hospitalar]
        destino = sorteio(pmf, espaco_amostral)
        preso.proximo_estado = state_name2code[destino]
        if destino == "recuperado"
            espaco_amostral = collect(simulacao.dias_graves_para_recuperados[1]:1:simulacao.dias_graves_para_recuperados[3])
            sdt1emeio = length(espaco_amostral) / 3
            pmf = exp.(-((espaco_amostral .- simulacao.dias_graves_para_recuperados[2]) .^ 2) / sdt1emeio)
            pmf = pmf / sum(pmf)
            dias_para_mudar = sorteio(pmf, espaco_amostral)
            preso.alterar_apos = dias_para_mudar
        elseif destino == "critico"
            espaco_amostral = collect(simulacao.dias_graves_para_criticos[1]:1:simulacao.dias_graves_para_criticos[3])
            sdt1emeio = length(espaco_amostral) / 3
            pmf = exp.(-((espaco_amostral .- simulacao.dias_graves_para_criticos[2]) .^ 2) / sdt1emeio)
            pmf = pmf / sum(pmf)
            preso.alterar_apos = sorteio(pmf, espaco_amostral)
        else
            espaco_amostral = collect(simulacao.dias_criticos_para_mortos[1]:1:simulacao.dias_criticos_para_mortos[3])
            sdt1emeio = length(espaco_amostral) / 3
            pmf = exp.(-((espaco_amostral .- simulacao.dias_criticos_para_mortos[2]) .^ 2) / sdt1emeio)
            pmf = pmf / sum(pmf)
            preso.alterar_apos = sorteio(pmf, espaco_amostral)
        end
        push!(simulacao.hospital.presos, preso)
        push!(simulacao.hospital.tempo_desde_chegada, 0)

        # Adicionar o custo de entrada no hospital
        simulacao.custo_total += custo_entrada_hospital
    end
    simulacao.prisao.presos = [simulacao.prisao.presos[i] for i in 1:length(simulacao.prisao.presos) if !(i in idx_agravar)]
end

function graves_para_criticos(simulacao, identificador)
    idx_criticar = findall(a -> (a.proximo_estado == state_name2code["critico"]) & (a.dia_alteracao_estado + a.alterar_apos == simulacao.dia), simulacao.hospital.presos)
    #open("entradas_na_uti_por_dia_$identificador.txt", "a") do file
        #write(file, string(length(idx_criticar)) * "\n")
    #end
    for i in idx_criticar
        preso = simulacao.hospital.presos[i]
        preso.estado = state_name2code["critico"]
        preso.dia_alteracao_estado = simulacao.dia
        espaco_amostral = ["recuperado", "morto_covid"]
        pmf = [1 - simulacao.mortalidade_hospitalar, simulacao.mortalidade_hospitalar]
        pmf = pmf .* (pmf .>= 0)
        pmf = pmf / sum(pmf)
        destino = sorteio(pmf, espaco_amostral)
        preso.proximo_estado = state_name2code[destino]
        if destino == "morto_covid"
            espaco_amostral = collect(simulacao.dias_criticos_para_mortos[1]:1:simulacao.dias_criticos_para_mortos[3])
            sdt1emeio = length(espaco_amostral) / 3
            pmf = exp.(-((espaco_amostral .- simulacao.dias_criticos_para_mortos[2]) .^ 2) / sdt1emeio)
            pmf = pmf / sum(pmf)
            dias_para_mudar = sorteio(pmf, espaco_amostral)
        elseif destino == "recuperado"
            espaco_amostral = collect(simulacao.dias_criticos_para_recuperados[1]:1:simulacao.dias_criticos_para_recuperados[3])
            sdt1emeio = length(espaco_amostral) / 3
            pmf = exp.(-((espaco_amostral .- simulacao.dias_criticos_para_recuperados[2]) .^ 2) / sdt1emeio)
            pmf = pmf / sum(pmf)
            dias_para_mudar = sorteio(pmf, espaco_amostral)
        end
        preso.alterar_apos = dias_para_mudar
    end
end

function morre_por_covid(simulacao, identificador)
    # Mortes no hospital
    idx_morre_hospital = findall(a -> (a.proximo_estado == state_name2code["morto_covid"]) && (a.dia_alteracao_estado + a.alterar_apos == simulacao.dia), simulacao.hospital.presos)
    for i in idx_morre_hospital
        preso = simulacao.hospital.presos[i]
        preso.dia_alteracao_estado = simulacao.dia
        preso.estado = state_name2code["morto_covid"]
        preso.alterar_apos = 0
        preso.proximo_estado = nothing
        push!(simulacao.mortos.presos, preso)
    end
    simulacao.hospital.presos = [simulacao.hospital.presos[i] for i in 1:length(simulacao.hospital.presos) if !(i in idx_morre_hospital)]

    # Total de mortes
    m = length(idx_morre_hospital)
    open("mortes_por_dia_$identificador.txt", "a") do file
        write(file, string(m) * "\n")
    end
end

function recupera(simulacao, identificador)
    n_recuperados = 0

    #recuperacao na prisao
    idx_recupera = findall(a->(a.proximo_estado == state_name2code["recuperado"]) & (a.dia_alteracao_estado+a.alterar_apos == simulacao.dia), simulacao.prisao.presos)
    n_recuperados += length(idx_recupera)
    for i in idx_recupera
        preso = simulacao.prisao.presos[i]
        preso.dia_alteracao_estado = simulacao.dia
        preso.estado = state_name2code["recuperado"]
        espaco_amostral = collect(simulacao.dias_recuperados_para_sucetiveis[1]:1:simulacao.dias_recuperados_para_sucetiveis[3])
        sdt1emeio = length(espaco_amostral)/3
        pmf = exp.(-((espaco_amostral.-simulacao.dias_recuperados_para_sucetiveis[2]).^2)/sdt1emeio)
        pmf = pmf/sum(pmf)
        preso.alterar_apos = sorteio(pmf, espaco_amostral)
        preso.proximo_estado = state_name2code["suscetivel"]
    end

    #recuperacao no hospital
    idx_recupera = findall(a->(a.proximo_estado == state_name2code["recuperado"]) & (a.dia_alteracao_estado+a.alterar_apos == simulacao.dia), simulacao.hospital.presos)
    n_recuperados += length(idx_recupera)
    #open("recuperados_por_dia_$identificador.txt", "a") do file
        #write(file, string(n_recuperados)*"\n")
    #end
    for i in idx_recupera
        preso = simulacao.hospital.presos[i]
        preso.estado = state_name2code["recuperado"]
        preso.dia_alteracao_estado = simulacao.dia
        espaco_amostral = collect(simulacao.dias_recuperados_para_sucetiveis[1]:1:simulacao.dias_recuperados_para_sucetiveis[3])
        sdt1emeio = length(espaco_amostral)/3
        pmf = exp.(-((espaco_amostral.-simulacao.dias_recuperados_para_sucetiveis[2]).^2)/sdt1emeio)
        pmf = pmf/sum(pmf)
        preso.alterar_apos = sorteio(pmf, espaco_amostral)
        preso.proximo_estado = state_name2code["suscetivel"]
        push!(simulacao.prisao.presos, preso)
    end

    simulacao.hospital.presos = [simulacao.hospital.presos[i] for i in 1:length(simulacao.hospital.presos) if !(i in idx_recupera)]
    n_recuperados = length(findall(a-> a.estado == state_name2code["recuperado"], simulacao.prisao.presos))
    #open("recuperados_$identificador.txt", "a") do file
        #write(file, string(n_recuperados)*"\n")
    #end
end

function transicao_sucetiveis(simulacao, identificador)
    idx_suscetivel = findall(a -> (a.proximo_estado == state_name2code["suscetivel"]) && (a.dia_alteracao_estado + a.alterar_apos == simulacao.dia), simulacao.prisao.presos)
    for i in idx_suscetivel
        preso = simulacao.prisao.presos[i]
        preso.estado = state_name2code["suscetivel"]
        preso.dia_alteracao_estado = simulacao.dia
        preso.alterar_apos = 0
        preso.proximo_estado = nothing
    end
end

function vacinacao(simulacao, identificador)
    if simulacao.dia in simulacao.dias_de_vacinacao
        n_vacinados = simulacao.vacinacao_por_dia
        idx_sucetiveis = findall(a-> (a.estado == state_name2code["suscetivel"] || a.estado == state_name2code["vacinado"]) && (a.dia_alteracao_estado+340) >= simulacao.dia, simulacao.prisao.presos)
        n_vacinados = min(n_vacinados, length(idx_sucetiveis))
        aux = sortperm(rand(length(idx_sucetiveis)))
        aux = aux[1:n_vacinados]
        vacinar = idx_sucetiveis[aux]

        for i in vacinar
            preso = simulacao.prisao.presos[i]
            preso.estado = state_name2code["vacinado"]
            preso.dia_alteracao_estado = simulacao.dia
            espaco_amostral = collect(simulacao.dias_vacinados_para_sucetiveis[1]:1:simulacao.dias_vacinados_para_sucetiveis[3])
            sdt1emeio = length(espaco_amostral)/3
            pmf = exp.(-((espaco_amostral.-simulacao.dias_vacinados_para_sucetiveis[2]).^2)/sdt1emeio)
            pmf = pmf/sum(pmf)
            preso.alterar_apos = sorteio(pmf, espaco_amostral)
            preso.proximo_estado = state_name2code["suscetivel"]

            # Adicionar o custo de entrada para vacinação
            simulacao.custo_total += custo_entrada_vacinacao
        end

        n_vacinados_hoje = length(vacinar)
        #open("vacinacao_por_dia_$identificador.txt", "a") do file
            #write(file, string(n_vacinados_hoje) * "\n")
        #end
    #else
        #open("vacinacao_por_dia_$identificador.txt", "a") do file
            #write(file, string(0) * "\n")
        #end
    end
end

function inserir_novos_presos(simulacao, identificador)
    tx_entrada = simulacao.dia <= 26 * 7 ? simulacao.tx_entrada_inicial : simulacao.tx_entrada
    num_vacinados = 0
    for i in 1:tx_entrada
        cela = rand(1:simulacao.n_celas)
        preso = instanciar_Preso(cela, simulacao.distribuicao_etaria, simulacao.n_grupos)
        if simulacao.dia > 366 && rand() < simulacao.p_vacinacao
            preso.estado = state_name2code["vacinado"]
            preso.proximo_estado = state_name2code["suscetivel"]
            preso.dia_alteracao_estado = simulacao.dia
            espaco_amostral = collect(simulacao.dias_vacinados_para_sucetiveis[1]:1:simulacao.dias_vacinados_para_sucetiveis[3])
            sdt1emeio = length(espaco_amostral) / 3
            pmf = exp.(-((espaco_amostral .- simulacao.dias_vacinados_para_sucetiveis[2]) .^ 2) / sdt1emeio)
            pmf = pmf / sum(pmf)
            preso.alterar_apos = sorteio(pmf, espaco_amostral)
            num_vacinados += 1
            
            # Adicionar o custo de entrada para vacinação
            simulacao.custo_total += custo_entrada_vacinacao
        end
        push!(simulacao.prisao.presos, preso)
    end
    #registrar_vacinados(num_vacinados, identificador)
end

function saida_prisao(simulacao, identificador)
    tx_saida = simulacao.dia <= 26 * 7 ? simulacao.tx_saida_inicial : simulacao.tx_saida
    estados_permitidos_saida = [
        state_name2code["suscetivel"],
        state_name2code["recuperado"],
        state_name2code["vacinado"],
        state_name2code["exposto"],
        state_name2code["infectado"]
    ]
    
    # Filtra os presos que podem sair
    presos_permitidos = [i for i in 1:length(simulacao.prisao.presos) if simulacao.prisao.presos[i].estado in estados_permitidos_saida]
    n_saem = min(tx_saida, length(presos_permitidos))  # Garante que não saia mais do que o número de presos permitidos
    
    if n_saem > 0
        aux = sortperm(rand(length(presos_permitidos)))
        saem = presos_permitidos[aux[1:n_saem]]
        
        # Atualiza o estado dos presos que estão saindo
        for idx in saem
            simulacao.prisao.presos[idx].estado = state_name2code["saida"]
        end

        # Registra o número de presos que saem
        #open("saida_$identificador.txt", "a") do file
            #write(file, string(n_saem) * "\n")
        #end
    end
end

# Função para calcular dias perdidos baseado no tempo em estados críticos
function calcular_dias_perdidos(preso)
    dias_perdidos = 0
    
    # Supondo que dias perdidos são os dias em que o preso esteve em estado grave ou crítico
    if preso.estado == state_name2code["grave"] || preso.estado == state_name2code["critico"] || preso.estado == state_name2code["morto_covid"]
        dias_perdidos += preso.dias_vividos  # Considera todos os dias vividos como dias perdidos
    end
    
    return dias_perdidos
end

# Função para calcular dias vividos com base na expectativa de vida e nos dias perdidos
function calcular_dias_vividos(preso)
    dias_perdidos = calcular_dias_perdidos(preso)
    dias_vividos = expectativa_vida_dias - dias_perdidos
    return dias_vividos
end

# Atualizar os dias vividos para cada preso na prisão e no hospital
function atualizar_dias_vividos(simulacao, identificador)
    for preso in simulacao.prisao.presos
        preso.dias_vividos = calcular_dias_vividos(preso)
    end

    for preso in simulacao.hospital.presos
        preso.dias_vividos = calcular_dias_vividos(preso)
    end
end

function contar_pessoas_na_prisao(simulacao, identificador)
    # Estados que devem ser incluídos na contagem
    estados_incluidos = [
        state_name2code["suscetivel"],
        state_name2code["exposto"],
        state_name2code["infectado"],
        state_name2code["vacinado"],
        state_name2code["recuperado"]
    ]
    
    # Contar apenas as pessoas que estão nos estados incluídos
    n_pessoas_prisao = count(p -> p.estado in estados_incluidos, simulacao.prisao.presos)
    
    return n_pessoas_prisao
end

function simular_dia(simulacao, identificador)
    # Inserir novos presos diariamente
    inserir_novos_presos(simulacao, identificador)
    
    # Saída de presos diariamente
    saida_prisao(simulacao, identificador)

    # Simulação das atividades diárias
    simular_celas(simulacao, identificador)
    simular_banho_de_sol(simulacao, identificador)

    # Transições de estado
    transicao_sucetiveis(simulacao, identificador)
    expostos_para_infectados(simulacao, identificador)
    infectados_para_graves(simulacao, identificador)
    graves_para_criticos(simulacao, identificador)
    morre_por_covid(simulacao, identificador)
    recupera(simulacao, identificador)

    # Vacinação (após todas as outras transições)
    vacinacao(simulacao, identificador)

    # Atualizar dias vividos
    atualizar_dias_vividos(simulacao, identificador)    

    # Registro dos estados (após todas as transições e vacinação)
    #registrar_suscetiveis(simulacao, identificador)
    #registrar_expostos(simulacao, identificador)
    #registrar_infectados(simulacao, identificador)
    
    # Calcular e registrar custo diário
    custo_diario = calcular_custo_diario(simulacao, identificador)
    simulacao.custo_total += custo_diario
    #registrar_custos_diarios(custo_diario, identificador)

    # Registro dos estados
    #n_pessoas_hospital = length(simulacao.hospital.presos)
    #open("pessoas_no_hospital_por_dia_$identificador.txt", "a") do file
        #write(file, string(n_pessoas_hospital) * "\n")
    #end

    #idx_uti = findall(a -> a.estado == state_name2code["critico"], simulacao.hospital.presos)
    #n_pessoas_uti = length(idx_uti)

    #open("pessoas_na_uti_por_dia_$identificador.txt", "a") do file
        #write(file, string(n_pessoas_uti) * "\n")
    #end

    #n_pessoas_prisao = contar_pessoas_na_prisao(simulacao, identificador)
    #open("pessoas_na_prisao_por_dia_$identificador.txt", "a") do file
        #write(file, string(n_pessoas_prisao) * "\n")
    #end

    # Registro dos vacinados após todas as mudanças do dia
    #n_pessoas_vacinado = length(findall(a -> a.estado == state_name2code["vacinado"], simulacao.prisao.presos))
    #open("vacinados_$identificador.txt", "a") do file
        #write(file, string(n_pessoas_vacinado) * "\n")
    #end

    simulacao.dia += 1
end
#=
function registrar_vacinados(num_vacinados, identificador)
    open("vacinados_entrando_por_dia_$identificador.txt", "a") do file
        write(file, string(num_vacinados) * "\n")
    end
end
=#
#=
function somatorio_dias_vividos(simulacao)
    # Calcula o somatório dos dias vividos por todos os presos
    total_dias_vividos = sum(preso -> preso.dias_vividos, simulacao.prisao.presos, init=0)
    total_dias_vividos += sum(preso -> preso.dias_vividos, simulacao.hospital.presos, init=0)
    total_dias_vividos += sum(preso -> preso.dias_vividos, simulacao.mortos.presos, init=0)
    return total_dias_vividos
end
=#
#=
function registrar_somatorio_dias_vividos(somatorio, identificador)
    open("somatorio_dias_vividos_$identificador.txt", "w") do file
        write(file, string(somatorio) * "\n")
    end
end
=#
function calcular_media_dias_vividos(simulacao)
    total_dias_vividos = sum(preso -> preso.dias_vividos, simulacao.prisao.presos, init=0)
    total_dias_vividos += sum(preso -> preso.dias_vividos, simulacao.hospital.presos, init=0)
    total_dias_vividos += sum(preso -> preso.dias_vividos, simulacao.mortos.presos, init=0)

    media_dias_vividos = total_dias_vividos / simulacao.total_pessoas
    return media_dias_vividos
end

#=
function registrar_media_dias_vividos(simulacao, identificador)
    media_dias_vividos = calcular_media_dias_vividos(simulacao)
    open("media_dias_vividos_$identificador.txt", "w") do file
        write(file, string(media_dias_vividos) * "\n")
    end
end
=#

function simular_semana(simulacao, identificador)
    for dia in 1:7
        simular_dia(simulacao, identificador)
    end
end
#=
# Função para registrar o AVG em um arquivo específico
function registrar_avg(avg, identificador, caminho_arquivo)
    open(caminho_arquivo, "a") do file
        write(file, string(identificador) * "\t" * string(avg) * "\n")
    end
end
=#
function simulacao_completa(parametros::Dict{String, Any}, identificador::String)
    # Extrair os parâmetros necessários
    tempo_de_simulacao = parametros["tempo_de_simulacao"]
    n_presos = parametros["n_presos"]
    tx_entrada = parametros["tx_entrada"] 
    tx_saida = parametros["tx_saida"]
    tx_entrada_inicial = parametros["tx_entrada_inicial"]
    tx_saida_inicial = parametros["tx_saida_inicial"]
    mortalidade_hospitalar = parametros["mortalidade_hospitalar"]
    n_celas = parametros["n_celas"]
    dias_de_vacinacao = collect(365:1:365+30)
    for i in 2:4
        dias_de_vacinacao = vcat(dias_de_vacinacao, collect(365*i:1:365*i+30))
    end
    tempo_sol = parametros["tempo_sol"]
    probabilidade_contagio = parametros["probabilidade_contagio"]
    tx_infectado_grave = parametros["tx_infectado_grave"]
    tx_grave_critico = parametros["tx_grave_critico"]
    vacinacao_por_dia = parametros["vacinacao_por_dia"]
    distribuicao_etaria = parametros["distribuicao_etaria"]
    distribuicao_etaria = distribuicao_etaria / sum(distribuicao_etaria)
    n_grupos = parametros["n_grupos"]
    tx_assintomaticos = parametros["tx_assintomaticos"]
    taxa_encontros_celas = parametros["taxa_encontros_celas"]
    p_vacinacao = parametros["p_vacinacao"]
    dias_recuperados_para_sucetiveis = parametros["dias_recuperados_para_sucetiveis"]
    dias_vacinados_para_sucetiveis = parametros["dias_vacinados_para_sucetiveis"]
    dias_expostos_para_infectados = parametros["dias_expostos_para_infectados"]
    dias_infectados_para_graves = parametros["dias_infectados_para_graves"]
    dias_graves_para_criticos = parametros["dias_graves_para_criticos"]
    dias_criticos_para_mortos = parametros["dias_criticos_para_mortos"]
    dias_infectados_para_recuperados = parametros["dias_infectados_para_recuperados"]
    dias_graves_para_recuperados = parametros["dias_graves_para_recuperados"]
    dias_criticos_para_recuperados = parametros["dias_criticos_para_recuperados"]
    total_pessoas = parametros["total_pessoas"]
    
    # Inicialização da simulação conforme seu código original
    prisao = instanciar_prisao(tempo_sol, n_celas, n_presos, distribuicao_etaria, n_grupos)
    hospital = Hospital([], [])
    mortos = Mortos([])
    custo_total = 0

    # Definir nomes de arquivos com base no identificador
    nomes_arquivos = [
        "casos_por_dia_$identificador.txt",
        "mortes_por_dia_$identificador.txt",
        "idade_dos_mortos_$identificador.csv"
    ]
    
    # Criar arquivos vazios ou diretórios conforme necessário
    for nome in nomes_arquivos
        open(nome, "w") do arquivo
        end
    end

    simulacao = Simulacao(
    tempo_de_simulacao,              # Int
    n_presos,                        # Int
    prisao,                          # Prisao
    hospital,                        # Hospital
    mortos,                          # Mortos
    tx_entrada,                      # Int
    tx_saida,                        # Int
    tx_entrada_inicial,              # Int
    tx_saida_inicial,                # Int
    mortalidade_hospitalar,          # Float64
    tx_infectado_grave,              # Vector{Float64}
    tx_grave_critico,                # Float64
    vacinacao_por_dia,               # Int
    1,                               # Int (dia inicial)
    probabilidade_contagio,          # Float64
    n_celas,                         # Int
    dias_de_vacinacao,               # Vector{Int}
    distribuicao_etaria,             # Vector{Float64}
    n_grupos,                        # Int
    tx_assintomaticos,               # Float64
    taxa_encontros_celas,            # Float64
    p_vacinacao,                     # Float64
    dias_recuperados_para_sucetiveis, # Vector{Float64}
    dias_vacinados_para_sucetiveis,  # Vector{Float64}
    dias_expostos_para_infectados,   # Vector{Float64}
    dias_infectados_para_graves,     # Vector{Float64}
    dias_graves_para_criticos,       # Vector{Float64}
    dias_criticos_para_mortos,       # Vector{Float64}
    dias_infectados_para_recuperados, # Vector{Float64}
    dias_graves_para_recuperados,    # Vector{Float64}
    dias_criticos_para_recuperados,  # Vector{Float64}
    custo_total,                     # Float64
    total_pessoas                    # Int
)


    simulacao.prisao.presos[1].estado = state_name2code["infectado"]
    simulacao.prisao.presos[1].proximo_estado = state_name2code["recuperado"]
    simulacao.prisao.presos[1].alterar_apos = 14

    for i in 1:simulacao.tempo_de_simulacao
        println("Semana:", i)
        simular_semana(simulacao, identificador)
        if i % 52 == 0
            for preso in simulacao.prisao.presos
                preso.idade += 1
            end

            for preso in simulacao.hospital.presos
                preso.idade += 1
            end
        end
    end

    # Calcular e registrar custo total
    #registrar_custo_total(simulacao.custo_total, identificador)

    # Aplicar taxa de desconto e registrar custo total descontado
    custo_total_descontado = aplicar_taxa_desconto(simulacao.custo_total, tempo_de_simulacao, 0.05)
    #registrar_custo_total_descontado(custo_total_descontado, identificador)

    # Calcular e registrar custo descontado por indivíduo
    custo_descontado_por_individuo = calcular_custo_descontado_por_individuo(custo_total_descontado, simulacao)
    #registrar_custo_descontado_por_individuo(custo_descontado_por_individuo, identificador)

    # Calcular média dos dias vividos
    media_dias_vividos = calcular_media_dias_vividos(simulacao)
    #registrar_media_dias_vividos(simulacao, identificador)

    # Calcular e registrar AVG
    avg = media_dias_vividos * 3 / 1092
    #registrar_avg(avg, identificador)

    # Calcular e registrar o somatório dos dias vividos
    #total_dias_vividos = somatorio_dias_vividos(simulacao)
    #registrar_somatorio_dias_vividos(total_dias_vividos, identificador)

    return simulacao
end

# Função para ler parâmetros de um arquivo e armazená-los em um dicionário
function ler_parametros(caminho_arquivo)
    variaveis = Dict{String, Any}()
    
    open(caminho_arquivo) do arquivo
        for linha in eachline(arquivo)
            if isempty(strip(linha)) || startswith(strip(linha), "#")
                continue
            end
            nome, valor_str = split(linha, '=')
            nome = strip(nome)
            valor_str = strip(valor_str)
            
            if startswith(valor_str, "[") && endswith(valor_str, "]")
                elementos = split(strip(valor_str, ['[', ']']), ',')
                valor = [parse(Float64, strip(elem)) for elem in elementos]
            else
                if occursin(".", valor_str)
                    valor = parse(Float64, valor_str)
                else
                    valor = parse(Int, valor_str)
                end
            end
            variaveis[nome] = valor
        end
    end
    
    return variaveis
end

# Lista de parâmetros e identificadores de cenário
cenarios = [
    ("cenario1", ler_parametros("parametros_cenario1.txt")),
    ("cenario2", ler_parametros("parametros_cenario2.txt")),
    ("cenario3", ler_parametros("parametros_cenario3.txt"))
]

# Defina um dicionário para armazenar os resultados dos três cenários
resultados = Dict{String, Dict{String, Float64}}()

# Caminho para os arquivos de custo-efetividade
caminho_custo_efetividade_1 = "custo_efetividade_1.txt"
caminho_custo_efetividade_2 = "custo_efetividade_2.txt"
caminho_avg_cenario1 = "AVG_cenario1.txt"
caminho_avg_cenario2 = "AVG_cenario2.txt"
caminho_avg_cenario3 = "AVG_cenario3.txt"
caminho_custo_descontado_cenario1 = "custo_descontado_individuo_cenario1.txt"
caminho_custo_descontado_cenario2 = "custo_descontado_individuo_cenario2.txt"
caminho_custo_descontado_cenario3 = "custo_descontado_individuo_cenario3.txt"

# Verificar se o arquivo já existe para evitar sobrescrever o cabeçalho
if !isfile(caminho_custo_efetividade_1)
    open(caminho_custo_efetividade_1, "w") do file
        write(file, "identificador\tdcusto1\tdavg1\n")
    end
end

if !isfile(caminho_custo_efetividade_2)
    open(caminho_custo_efetividade_2, "w") do file
        write(file, "identificador\tdcusto2\tdavg2\n")
    end
end

# Função para registrar o AVG em um arquivo específico
function registrar_avg(avg, identificador, caminho_arquivo)
    open(caminho_arquivo, "a") do file
        write(file, string(identificador) * "\t" * string(avg) * "\n")
    end
end

# Função para registrar o custo descontado por indivíduo em um arquivo específico
function registrar_custo_descontado(custo_descontado_individuo, identificador, caminho_arquivo)
    open(caminho_arquivo, "a") do file
        write(file, string(identificador) * "\t" * string(custo_descontado_individuo) * "\n")
    end
end
# Executar simulação para cada cenário
for (nome_cenario, variaveis) in cenarios
    println("Iniciando simulação para $nome_cenario")

    # Gerar um identificador único combinando o timestamp e o nome do cenário
    identificador = string(nome_cenario, "_", time)

    # Executar a simulação
    simulacao_teste = simulacao_completa(variaveis, identificador)

    # Calcular custo descontado por indivíduo e avg
    custo_total_descontado = aplicar_taxa_desconto(simulacao_teste.custo_total, variaveis["tempo_de_simulacao"], 0.05)
    custo_descontado_por_individuo = calcular_custo_descontado_por_individuo(custo_total_descontado, simulacao_teste)
    media_dias_vividos = calcular_media_dias_vividos(simulacao_teste)
    avg = media_dias_vividos * 3 / 1092

    # Armazenar os resultados para cada cenário
    resultados[nome_cenario] = Dict(
        "custo_descontado_por_individuo" => custo_descontado_por_individuo,
        "avg" => avg
    )

    # Criar DataFrame para armazenar as idades dos mortos
    df = DataFrame(idade=Int[])
    for morto in simulacao_teste.mortos.presos
        push!(df, (idade = morto.idade,))
    end

    # Salvar o DataFrame em um arquivo CSV
    CSV.write("idade_dos_mortos_$identificador.csv", df)

    println("Simulação para $nome_cenario concluída e arquivo idade_dos_mortos_$identificador.csv gerado.")

    # Determinar o caminho do arquivo AVG com base no cenário
    caminho_avg = if nome_cenario == "cenario1"
        caminho_avg_cenario1
    elseif nome_cenario == "cenario2"
        caminho_avg_cenario2
    else
        caminho_avg_cenario3
    end

    # Registrar o AVG no arquivo correspondente
    registrar_avg(avg, identificador, caminho_avg)
    # Calcular custo descontado por indivíduo
    custo_total_descontado = aplicar_taxa_desconto(simulacao_teste.custo_total, variaveis["tempo_de_simulacao"], 0.05)
    custo_descontado_por_individuo = calcular_custo_descontado_por_individuo(custo_total_descontado, simulacao_teste)

    # Determinar o caminho do arquivo de custo descontado com base no cenário
    caminho_custo_descontado = if nome_cenario == "cenario1"
        caminho_custo_descontado_cenario1
    elseif nome_cenario == "cenario2"
        caminho_custo_descontado_cenario2
    else
        caminho_custo_descontado_cenario3
    end

    # Registrar o custo descontado no arquivo correspondente
    registrar_custo_descontado(custo_descontado_por_individuo, identificador, caminho_custo_descontado)

    # Após a execução de cada cenário, calcular e salvar os arquivos de custo-efetividade
    if all(key -> haskey(resultados, key), ["cenario1", "cenario2", "cenario3"])
        # Calcular dcusto1, davg1, dcusto2, davg2
        dcusto1 = resultados["cenario2"]["custo_descontado_por_individuo"] - resultados["cenario1"]["custo_descontado_por_individuo"]
        davg1 = resultados["cenario2"]["avg"] - resultados["cenario1"]["avg"]

        dcusto2 = resultados["cenario3"]["custo_descontado_por_individuo"] - resultados["cenario1"]["custo_descontado_por_individuo"]
        davg2 = resultados["cenario3"]["avg"] - resultados["cenario1"]["avg"]

        identificador1 = time
        
        # Adicionar os resultados ao arquivo custo_efetividade_1
        open(caminho_custo_efetividade_1, "a") do file
            write(file, string(identificador1) * "\t" * string(dcusto1) * "\t" * string(davg1) * "\n")
        end

        # Adicionar os resultados ao arquivo custo_efetividade_2
        open(caminho_custo_efetividade_2, "a") do file
            write(file, string(identificador1) * "\t" * string(dcusto2) * "\t" * string(davg2) * "\n")
        end

        println("Resultados de custo-efetividade adicionados aos arquivos custo_efetividade_1.txt e custo_efetividade_2.txt.")
    else
        println("Erro: Certifique-se de que os três cenários foram simulados corretamente.")
    end
end

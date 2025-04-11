__precompile__() # Este comando es para que julia precompile el paquete


module quantum

using LinearAlgebra 

export random_state, projector, apply_unitary
export apply_unitary!, partial_trace
export I, X, Y, Z, sigmas, sigma

@doc "Función para connvertir una lista de bits a un número entero."
function original_integer(list)
    return parse(Int, join(list); base=2)
end

@doc "Función para convertir un número entero a una lista de bits. 
    Si se especifica el argumento pad, la lista de bits tendrá longitud pad."
function base_2(integer; pad= nothing)
    if pad == nothing
        return reverse(digits(integer, base = 2))
    else
        return reverse(digits(integer, base = 2, pad = pad))
    end
end

@doc "La mascara toma un número entero que indica a cuales qubits se les aplicará la operación. 
Por ejemplo, si mask=4, indicará que se le aplicará la operación al tercer qubit [1,0,0].

a es el targert (sub espacio de qubits a los que se les aplica la operación)
y b el untouched  (sub espacio de qubits que no se les aplica la operación).

La función regresa el indice de la compenente del estado que se le aplica la operación."
function merge_two_integers(mask, a, b, n)
    result = 0
    i=0
    bit_mask= base_2(mask, pad=n)
    while i<=n-1 
        if bit_mask[end-i] == 1
            bit = a & 1
            result = result | bit<<i
            a = a>>1
        else
            bit = b & 1
            result = result | bit<<i
            b = b>>1
        end
        i+=1
    end

    return result
               
end


@doc "Función para generar un estado cuántico aleatorio de dimensión dim. 
    El estado es un vector columna normalizado."
function random_state(dim::Int=2)::Vector{Complex{<:AbstractFloat}}
    v=randn(dim,1)+randn(dim,1)im
    v=v/norm(v)
    return vec(v)
end

@doc "Función para generar el proyector de un estado cuántico. 
    El proyector es una matriz que representa la proyección del estado en sí mismo."
function projector(state::Vector{T}) where T
    return state*state'
end

@doc "Función para aplicar una operación unitaria a un estado cuántico. 

    U es la operación unitaria que se aplica en solo un subespacio del sistema completo.

    La operación se aplica a un subespacio de qubits especificado por el argumento target. Si target=3, 
    [0, 1, 1] indica que se aplica al 1 y 2 qubit.


    La función regresa el nuevo estado cuántico."
function apply_unitary(state::Vector{T}, U, target) where T
    n=Int(log2(length(state)))
    new_state = zeros(ComplexF64, size(state))

    bit_mask = base_2(target)
    n_t= count(x->x==1, bit_mask)
    target_dim= 2^n_t
    unt_dim= Int(length(state)/target_dim)
    for j in 0:unt_dim-1
        pos=[merge_two_integers(target, i, j,n) + 1 for i in 0:target_dim-1]
        new_state[pos]=U*state[pos]
    end
    return new_state
end

function apply_unitary(state::Matrix{T}, U, target) where T
    new_state = zeros(ComplexF64, size(state))
    n=Int(log2(size(state, 1)))

    bit_mask= base_2(target)
    n_t=count(x->x==1, bit_mask)
    
    target_dim= 2^n_t
    unt_dim_j= Int(size(state, 1)/target_dim)
    unt_dim_k= Int(size(state, 2)/target_dim)
    
    for j in 0:unt_dim_j-1
        for k in 0:unt_dim_k-1
            pos_j=[merge_two_integers(target, i, j, n) + 1 for i in 0:target_dim-1]
            pos_k=[merge_two_integers(target, i, k, n) + 1 for i in 0:target_dim-1]
            new_state[pos_j, pos_k]= U*state[pos_j, pos_k]*U'
        end
    end
    return new_state
end

function apply_unitary!(state::Vector{T}, U, target) where T
    n=Int(log2(length(state)))
    
    bit_mask = base_2(target)
    n_t= count(x->x==1, bit_mask)
    target_dim= 2^n_t
    unt_dim= Int(length(state)/target_dim)
    for j in 0:unt_dim-1
        pos=[merge_two_integers(target, i, j,n) + 1 for i in 0:target_dim-1]
        state[pos]=U*state[pos]
    end
    
end

function apply_unitary!(state::Matrix{T}, U, target) where T
    
    n=Int(log2(size(state, 1)))

    bit_mask= base_2(target)
    n_t=count(x->x==1, bit_mask)
    
    target_dim= 2^n_t
    unt_dim_j= Int(size(rho, 1)/target_dim)
    unt_dim_k= Int(size(rho, 2)/target_dim)
    
    for j in 0:unt_dim_j-1
        for k in 0:unt_dim_k-1
            pos_j=[merge_two_integers(target, i, j, n) + 1 for i in 0:target_dim-1]
            pos_k=[merge_two_integers(target, i, k, n) + 1 for i in 0:target_dim-1]
            state[pos_j, pos_k]= U*state[pos_j, pos_k]*U'
        end
    end
    
end

@doc "Función para calcular la traza parcial de un estado cuántico. 

    El argumento target indica cuáles qubits no van a ser trazados. Si target= 1, [0,0, 1] indica que 
    se va a trazar sobre el subespacio del segundo y tercer qubit."
function partial_trace(state::Matrix{T}, target) where T
    n=Int(log2(size(state, 1)))
    bit_mask = base_2(target)
    n_t= count(x->x==1, bit_mask)

    target_dim=2^n_t
    n_unt= n-n_t
    unt_dim= 2^n_unt

    new_state = zeros(ComplexF64, target_dim, target_dim)
    for  i in 0:target_dim-1
        for j in 0:target_dim-1
            for k in 0:unt_dim-1
                pos_i=merge_two_integers(target, i, k, n) + 1 
                pos_j=merge_two_integers(target, j, k, n) + 1 
                new_state[i+1, j+1] += state[pos_i, pos_j]

            end
        end
    end
    return new_state
end


I=[ 1.0 0.0; 0.0 1.0]
X=[0.0 1.0; 1.0 0.0]
Y=[0.0 -1.0im; 1.0im 0.0]
Z=[1.0 0.0; 0.0 -1.0]

@doc "Diccionario de matrices sigma. 
    1: Pauli-X, 2: Pauli-Y, 3: Pauli-Z, 0: Identidad."
sigmas=Dict(1=>X, 2=>Y, 3=>Z, 0=>I)


@doc "Indice, indica la matriz de Pauli que se va a usar. 
    pos, indica la posición de la matriz en el tensor producto. 
    n, indica el número de matrices en el tensor producto."
function sigma(indice, pos, n)
    mat=sigmas[indice]
    list=[]
    for i in 1:n
        if i==pos
            push!(list, mat)
        else
            push!(list, sigmas[0])
        end
    end
    return kron(list ...)
end

function Ising(n, b, J; cerrada=false)
    if cerrada==false
        m=n-1
       
    elseif cerrada==true
        m=n
        
    end
   

    H=-sum([J*sigma(3, mod(i, n+1), n)*sigma(3,mod(i+1, n+1), n) for i in 1:m]) + b*sum([sigma(1, i, n) for i in 1:n])

    
    return H
    
end


end
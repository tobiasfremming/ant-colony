
module MyDataStructures
    using JSON3

    struct Patient
        id::Int
        demand::Float64
        start_time::Float64
        end_time::Float64
        care_time::Float64
        x_coord::Float64
        y_coord::Float64
    end

    mutable struct Nurse
        id::Int
        capacity::Float64
        route::Vector{Int}  # Represents a sequence of patients (permutation)
    end

    mutable struct Instance
        nurses::Vector{Nurse}
        patients::Vector{Patient}
        travel_times::Matrix{Float64}
        depot_x::Float64
        depot_y::Float64
        depot_return_time::Float64
        fitness::Float64
    end

    mutable struct Population
        instances::Vector{Instance}  # List of different solutions (GA population)
    end


    function parseInstance(file_path::String)::Instance
        if !isfile(file_path)
            error("File not found: $file_path. Ensure the JSON file is in the correct directory.")
        end
        json_data = JSON3.read(read(file_path, String))

        # Convert patient IDs from symbols to integers
        patients = [Patient(parse(Int, String(k)), v[:demand], v[:start_time], v[:end_time], v[:care_time], v[:x_coord], v[:y_coord]) 
                    for (k, v) in pairs(json_data[:patients])]

        nurses = [Nurse(i, json_data[:capacity_nurse], []) for i in 1:json_data[:nbr_nurses]]

        # Convert JSON3 nested array to a standard Julia matrix
        travel_times = reduce(hcat, [collect(Float64, row) for row in json_data[:travel_times]])'

        depot_x = json_data[:depot][:x_coord]
        depot_y = json_data[:depot][:y_coord]
        depot_return_time = json_data[:depot][:return_time]

        return Instance(nurses, patients, travel_times, depot_x, depot_y, depot_return_time, 1.)
    end





end


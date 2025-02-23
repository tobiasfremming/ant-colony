function read_from_json(file_name)
    open(file_name) do file
        return JSON.parse(file)
    end
    
end

function write_to_json(file_name, data)
    open(file_name, "w") do file
        JSON.print(file, data)
    end
end



# function calculate_travel_time(ant::Ant, travel_time_matrix::Matrix{Float64}, depot::Depot, food_sources::Vector{FoodSource})
#     total_travel_time = 0.0
#     current_node = 1  # Depot index
#     current_time = 0.0 

#     for food_source_index in ant.chromosome
#         food_source = food_sources[food_source_index]  # Get patient data
#         travel_time = travel_time_matrix[current_node, food_source_index + 1]  # +1 because index 1 is depot
#         total_travel_time += travel_time
#         current_time += travel_time
#         current_time += food_source.care_time
#         current_x, current_y = food_source.x_coord, food_source.y_coord
#     end

#     travel_time = travel_time_matrix[current_node, 1]  # Return to depot
#     total_travel_time += travel_time
#     current_time += travel_time

#     return total_travel_time
# end

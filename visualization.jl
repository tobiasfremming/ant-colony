module Visualization
    import ..MyDataStructures: Instance, Patient, Nurse 
    using Plots

    function plotNurseRoutes(instance::Instance)
        plt = plot()
        
        # Plot depot at correct location
        scatter!([instance.depot_x], [instance.depot_y], markersize=8, label="Depot", color=:black)
        
        for nurse in instance.nurses
            if isempty(nurse.route)
                continue
            end
            
            x_coords = [instance.depot_x]  # Start at depot
            y_coords = [instance.depot_y]
            
            for patient_id in nurse.route
                patient = instance.patients[patient_id]
                push!(x_coords, patient.x_coord)
                push!(y_coords, patient.y_coord)
            end
            
            push!(x_coords, instance.depot_x)  # Return to depot
            push!(y_coords, instance.depot_y)
            
            plot!(x_coords, y_coords, markershape=:circle, label="Nurse $(nurse.id)")
        end
        display(plt)
    end
end

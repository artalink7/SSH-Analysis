function SaveWithTimeStamp(data, base_dir::String, prefix::String)
    # Ensure the directory exists
    isdir(base_dir) || mkpath(base_dir)

    filename = joinpath(base_dir, "$(prefix).txt")
    writedlm(filename, data)
    println("Saved to: $filename")
end

function ProjectRoot()
    return normpath(joinpath(@__DIR__, ".."))
end

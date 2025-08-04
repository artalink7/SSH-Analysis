function SaveWithTimeStamp(data, base_dir::String, prefix::String)
    filename = joinpath(base_dir, "$(prefix).txt")
    writedlm(filename, data)
    println("Saved to: $filename")
end

function ProjectRoot()
    return normpath(joinpath(@__DIR__, ".."))
end
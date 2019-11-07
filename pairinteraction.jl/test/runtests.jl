detected_tests = filter(
    name->startswith(name, "test_") && endswith(name, ".jl"),
    readdir("."))

for name in detected_tests
    if startswith(name, "test_") && endswith(name, ".jl")
        include(name)
    end
end

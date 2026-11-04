// Resolve params.gurobi_license to a file, or [] when there is none to stage.
// Unset/blank -> []; the default name missing -> [] plus a warning (gurobipy uses its
// own lookup); any other missing path -> error.
def gurobiLicenseFile(value) {
    if (!value || value instanceof Boolean || value.toString() == "true") {
        return []
    }
    def lic = file(value.toString())
    if (lic.exists()) {
        return lic
    }
    if (value.toString() == "gurobi.lic") {
        log.warn "No gurobi.lic in the launch directory -- Gurobi (if used) falls back to its own license lookup. Pass --gurobi_license to stage one."
        return []
    }
    error("--gurobi_license '${value}' does not exist.")
}

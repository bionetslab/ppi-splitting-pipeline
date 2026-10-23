// Seeds for the three independent sources of randomness

def robustnessSeeds() {
    [split: [1234, 1, 2], sampling: [42, 1, 2], training: [42, 1, 2]]
}

def seedsFor(String level) {
    if (params.robustness_test) {
        return robustnessSeeds()[level]
    }
    return [[split: params.split_seed, sampling: params.sampling_seed, training: params.training_seed][level]]
}

// Fans every (meta, ...) item out over the seeds of the given levels
// ("split", "sampling", "training"), adding meta.<level>_seed. Expanding a
// channel to the same levels always yields the same metas, so channels
// expanded separately still join on meta. Outside --robustness_test every
// level has a single seed, so this is a 1:1 map.
def expandSeeds(ch, List levels) {
    ch.flatMap { item ->
        def metas = [item[0]]
        levels.each { level ->
            def key = level + '_seed'
            metas = metas.collectMany { m -> seedsFor(level).collect { s -> m + [(key): s] } }
        }
        metas.collect { m -> [m] + item.drop(1) }
    }
}

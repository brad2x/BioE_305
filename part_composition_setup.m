function part_composition_setup(version)
    main_url = 'http://web.mit.edu/20.305/www';
    files = {'BioSystem' 'Compositor' 'Const' 'Part' 'Pulse' 'Rate'};
    good_versions = {'v1pre' 'v2pre' 'v3' 'v4' 'v5' 'v6'};
    assert(ismember(version, good_versions), sprintf('Invalid version %s!', version))

    for f = files
        file = f{1};
        url = sprintf('%s/%s/%s.m', main_url, version, file);
        urlwrite(url, sprintf('%s.m', file));
    end

    rehash
end

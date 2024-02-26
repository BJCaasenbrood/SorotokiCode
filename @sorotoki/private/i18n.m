function str = i18n(key, varargin)
    persistent locale nls
    if ~isstruct(nls)
        locale = char(regexp(get(0, 'Language'), '^[a-zA-Z]+', 'match'));
        nls = load('mpi_nls.mat');
    end

    %% Check if message key exists.
    if isfield(nls, locale)
        data = nls.(locale);
    else
        data = nls.en;
    end

    if ~ischar(key)
        if (                                                                ...
            ~isfield(nls, locale)                                           ...
            || ~isfield(nls.(locale), 'unexpected_key')                     ...
        )
            error(nls.en.unexpected_key, class(key));
        else
            error(nls.(locale).unexpected_key, class(key));
        end
    end

    if (                                                                    ...
        ~isfield(data, key)                                                 ...
        && strcmp(locale, 'en')                                             ...
        || ~isfield(nls.en, key)                                            ...
    )
        error(nls.en.undefined_key, key);
    end

    %% Get the localised message.
    if isfield(data, key)
        str = data.(key);
    else
        str = sprintf(nls.en.(key), string(varargin{:}));
    end

    %% Variable argument substitution.
    if nargin == 1
        return;
    end

    xs = cellfun(@char, varargin, 'uni', false);
    str = sprintf(str, xs{:});
end

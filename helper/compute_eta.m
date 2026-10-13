function eta = compute_eta(n, S, constants)

    eta =[];
    for ss=1:S
        for kk=1:n
            eta(ss,kk) = sum(constants{ss}(:,:,kk));
        end
    end

end
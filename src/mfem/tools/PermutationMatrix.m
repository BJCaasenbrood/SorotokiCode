function PMat = PermutationMatrix(NElem)
if ispc
    usermem = memory; MaxArray = usermem.MaxPossibleArrayBytes / 8;
    if (NElem*NElem > MaxArray )
        cprintf('error','* not enough dynamic memory available!');
        cprintf('error','reduce the number of elements');
        PMat = 0;
    else
        PMat = VChooseK(1:NElem,2);
    end
else
    PMat = VChooseK(1:NElem,2);
end

end
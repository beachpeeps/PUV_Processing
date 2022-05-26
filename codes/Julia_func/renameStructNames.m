function newS = renameStructNames(oldS, oldNames, newNames)
newS = struct;
for whichField = 1:numel(oldNames)
  oldN = oldNames{whichField};
  if ~isfield(oldS, oldN)
    error(sprintf('Field %s does not exist in old struct array!', oldN));
  end
  newN = newNames{whichField};
  newS.(newN) = oldS.(oldN);
end
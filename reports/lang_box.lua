local labels = {
  r = "R",
  bash = "Bash",
  sh = "Bash",
  python = "Python"
}

function CodeBlock(el)
  local lang = el.classes[1]
  if lang == nil then return nil end
  local label = labels[string.lower(lang)]
  if label == nil then return nil end
  return {
    pandoc.RawBlock("latex", "\\begin{langbox}{" .. label .. "}"),
    el,
    pandoc.RawBlock("latex", "\\end{langbox}")
  }
end

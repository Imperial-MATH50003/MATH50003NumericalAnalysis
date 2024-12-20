\documentclass[12pt,a4paper]{article}

\usepackage[a4paper,text={16.5cm,25.2cm},centering]{geometry}
\usepackage{lmodern}
\usepackage{amssymb,amsmath}
\usepackage{bm}
\usepackage{graphicx}
\usepackage{microtype}
\usepackage{hyperref}
\usepackage[usenames,dvipsnames]{xcolor}
{{#:tex_deps}}
{{{ :tex_deps }}}
{{/:tex_deps}}
\setlength{\parindent}{0pt}
\setlength{\parskip}{1.2ex}




\hypersetup
       {   pdfauthor = { {{{:author}}} },
           pdftitle={ {{{:title}}} },
           colorlinks=TRUE,
           linkcolor=black,
           citecolor=blue,
           urlcolor=blue
       }

{{#:title}}
\title{ {{{ :title }}} }
{{/:title}}

{{#:author}}
\author{ {{{ :author }}} }
{{/:author}}

{{#:date}}
\date{ {{{ :date }}} }
{{/:date}}

{{ :highlight }}

\def\endash{–}
\def\bbD{ {\mathbb D} }
\def\bbZ{ {\mathbb Z} }
\def\bbR{ {\mathbb R} }
\def\bbC{ {\mathbb C} }

\def\x{ {\vc x} }
\def\a{ {\vc a} }
\def\b{ {\vc b} }
\def\e{ {\vc e} }
\def\f{ {\vc f} }
\def\u{ {\vc u} }
\def\v{ {\vc v} }
\def\y{ {\vc y} }
\def\z{ {\vc z} }
\def\w{ {\vc w} }

\def\bt{ {\tilde b} }
\def\ct{ {\tilde c} }
\def\Ut{ {\tilde U} }
\def\Qt{ {\tilde Q} }
\def\Rt{ {\tilde R} }
\def\Xt{ {\tilde X} }
\def\acos{ {\rm acos}\, }

\def\red#1{ {\color{red} #1} }
\def\blue#1{ {\color{blue} #1} }
\def\green#1{ {\color{ForestGreen} #1} }
\def\magenta#1{ {\color{magenta} #1} }

\input{somacros}

\begin{document}

{{#:title}}\maketitle{{/:title}}

{{{ :body }}}

\end{document}
% Richard Glennie latex document class file 
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{rgclass}[2019/04/30 Richard Glennie document class]

%% class inhertance
\LoadClass[twocolumn, twoside, a4paper, final]{article} 

%% options 
% referee option 
\newif\ifSFB@referee
\def\ds@referee{%
 \SFB@refereetrue
 \@twocolumnfalse
}
% pass all unknown options to article
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax

%% packages
% citations
\RequirePackage{natbib} 
% maths symbols
\RequirePackage{amsmath, amsthm, amssymb}
% graphics 
\RequirePackage{graphicx}
% hyperlinks 
\RequirePackage{hyperref}
% bold math
\RequirePackage{bm}
% url command 
\RequirePackage{url}
% colours 
\RequirePackage{xcolor}
% headers and footers 
\RequirePackage{fancyhdr}
% microtype: typographical fixes
\RequirePackage{microtype} 
% ensure correct hyphenation 
\RequirePackage[UKenglish]{babel}
% better section headers 
\RequirePackage[compact, small, sc]{titlesec}
% tables 
\RequirePackage{booktabs} 

%% Lines 
\lineskip 1pt
\normallineskip 1pt
\ifSFB@referee
  \def\baselinestretch{1.8}
\else
  \def\baselinestretch{1}
\fi

%% Page layout 
% 2 column by default, 1 column in referee mode 
\if@twocolumn
  \setlength\oddsidemargin{-2.5pc}
  \setlength\evensidemargin{-2.5pc}
\else
  \setlength\oddsidemargin{-0.15in}%\z@}
  \setlength\evensidemargin{-0.15in}
\fi
\setlength\marginparwidth{0.66in}

\ifSFB@referee
  \setlength\topmargin{-5.1pc}%
\else
  \setlength\topmargin{-2.5pc}%\z@}
\fi
\setlength{\headheight}{0.2in}
\setlength{\headsep}   {0.3in}
\setlength{\topskip}   {0.15in}
\setlength{\footskip}  {0.3in}

\setlength\textheight{9in}
\setlength\textwidth{7in}
% use fancy headers, empty by default 
\pagestyle{fancy}
\rhead{}
\lhead{}
\chead{}
% remove horizontal line below header 
\renewcommand{\headrulewidth}{0pt}

\ifSFB@referee
 \setlength\textwidth{7in}  
\fi
\setlength\columnsep{0.15in}
\setlength\columnseprule{0in}
\def\titlepage{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
  \else \newpage \fi \thispagestyle{empty}\c@page\z@}
\def\endtitlepage{\if@restonecol\twocolumn \else \newpage \fi}

%% Title page
% title page: title, authors, affiliations, date, abstract 
\RequirePackage{titling} 
\pretitle{\vspace*{23pt}\raggedright\sloppy\LARGE\begin{center}\bf}
\posttitle{\end{center}\vspace{1.5pt}} 
\preauthor{\large\begin{center}\lineskip .75em%
\begin{tabular}[t]{c}}
\postauthor{\end{tabular}\par\end{center}\vspace{1.5pt}}
\predate{\large\begin{center}}
\postdate{\end{center}\vspace{5pt}}
\thanksmarkseries{arabic}
\RequirePackage{authblk}
\setlength{\affilsep}{0.5em}
% abstract 
\renewcommand{\abstract}[2][.8]{
  \renewcommand\maketitlehookd{
    \mbox{}\medskip\par
    \centering
    \begin{minipage}{#1\textwidth}
        {\sc Summary:} #2
    \end{minipage}
    \vspace{18pt}
  }
} 


%% Fonts 
\if@twocolumn
  \renewcommand\normalsize{\@setfontsize\normalsize\@ixpt{11}%
   \abovedisplayskip 6pt plus 2pt minus 2pt
   \belowdisplayskip \abovedisplayskip
   \abovedisplayshortskip 6pt plus 2pt
   \belowdisplayshortskip \abovedisplayshortskip
   \let\@listi\@listI}
 \else
 \ifSFB@referee
  \renewcommand\normalsize{\@setfontsize\normalsize\@xiipt{14.5}% From size12.clo
   \abovedisplayskip 4pt plus 1pt minus 1pt
   \belowdisplayskip \abovedisplayskip
   \abovedisplayshortskip 4pt plus 1pt
   \belowdisplayshortskip \abovedisplayshortskip
   \let\@listi\@listI}
 \else
 \renewcommand\normalsize{\@setfontsize\normalsize\@ixpt{11}%
   \abovedisplayskip 4pt plus 1pt minus 1pt
   \belowdisplayskip \abovedisplayskip
   \abovedisplayshortskip 4pt plus 1pt
   \belowdisplayshortskip \abovedisplayshortskip
   \let\@listi\@listI}
 \fi
\fi
% Garamond by default 
\RequirePackage{eulervm}
\RequirePackage{ebgaramond}
\RequirePackage[T1]{fontenc}

%% Floats
\setlength\floatsep{0.16in}
\setlength\textfloatsep{0.16in}
\setlength\intextsep{0.16in}
\setlength\dblfloatsep{0.16in}
\setlength\dbltextfloatsep{0.15in}
\setlength\marginparpush{0.05in}

%% Paragraphs 
\newlength{\realparindent}
\setlength\parskip{ \z@}
\setlength\parindent{10pt}
\setlength{\realparindent}{\parindent}
\setlength\partopsep{0pt}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\@beginparpenalty -\@lowpenalty
\@endparpenalty -\@lowpenalty
\@itempenalty -\@lowpenalty
% penalities stop hanging lines 
\hyphenpenalty=500
\clubpenalty=10000
\widowpenalty=10000

%% Sections
\newcommand*{\Section}{}
\let\Section\section
\renewcommand*{\section}{%
  \@ifstar {\star@section}{\@dblarg\nonstar@section}%
}
\newcommand*{\star@section}[1]{%
  \Section*{\MakeLowercase{#1}}%
}
\newcommand*{\nonstar@section}[2][]{%
  \Section[{#1}]{\MakeLowercase{#2}}%
}
\let\Subsection\subsection
\renewcommand*{\subsection}{%
  \@ifstar {\star@subsection}{\@dblarg\nonstar@subsection}%
}
\newcommand*{\star@subsection}[1]{%
  \Subsection*{\MakeLowercase{#1}}%
}
\newcommand*{\nonstar@subsection}[2][]{%
  \Subsection[{#1}]{\MakeLowercase{#2}}%
}
\let\Subsubsection\subsubsection
\renewcommand*{\subsubsection}{%
  \@ifstar {\star@subsubsection}{\@dblarg\nonstar@subsubsection}%
}
\newcommand*{\star@subsubsection}[1]{%
  \Subsubsection*{\MakeLowercase{#1}}%
}
\newcommand*{\nonstar@subsubsection}[2][]{%
  \Subsubsection[{#1}]{\MakeLowercase{#2}}%
}
\def\thesection {\arabic{section}}
\def\thesubsection    {\thesection.\arabic{subsection}}
\def\thesubsubsection {\thesubsection.\arabic{subsubsection}}
\def\theparagraph     {\thesubsubsection.\arabic{paragraph}}
\def\thesubparagraph  {\theparagraph.\arabic{subparagraph}}

%% Lists
\setlength\leftmargini{   1em }
\setlength\leftmarginii{  1em }
\setlength\leftmarginiii{ 1em }
\setlength\leftmarginiv{  1em }
\setlength\leftmarginv{  1pc }
\setlength\leftmarginvi{  1pc }
\setlength\leftmargin{\leftmargini}
\setlength\labelsep{ 0.5em}
\labelwidth\leftmargini\advance\labelwidth-\labelsep
\def\labelitemi{$\bullet$}
\def\labelitemii{--}
\def\labelitemiii{$\cdot$}
\def\labelitemiv{*}

\def\labelenumi{\theenumi}
\def\theenumi{(\arabic{enumi})}% \roman
\def\labelenumii{\theenumii}
\def\theenumii{(\alph{enumii})}
\def\p@enumii{\theenumi}
\def\labelenumiii{\theenumiii}
\def\theenumiii{(\arabic{enumiii})}
\def\p@enumiii{\theenumi(\theenumii)}
\def\labelenumiv{\theenumiv}
\def\theenumiv{(\Alph{enumiv})}
\def\p@enumiv{\p@enumiii\theenumiii}

%% Quotations
\RequirePackage{csquotes}
\renewcommand{\mkbegdispquote}[2]{\itshape}

%% Custom environments
\theoremstyle{definition}
\newtheorem{defn}{Definition}[section] 
\theoremstyle{plain}
\newtheorem{thm}{Theorem}[section]
\newtheorem{prop}{Proposition}[section]
\newtheorem{coro}{Corollary}[section]
\newtheorem{lemma}{Lemma}[section]

%% Custom comamnds
\newcommand{\rmn}[1]{{\mathrm{#1}}}
\newcommand{\itl}[1]{{\mathit{#1}}}
\newcommand{\bld}[1]{{\mathbf{#1}}}
\newcommand{\be}{\begin{itemize}}
\newcommand{\en}{\end{itemize}}
\newcommand{\beN}{\begin{enumerate}}
\newcommand{\enN}{\end{enumerate}}

% end of rgclass.cls




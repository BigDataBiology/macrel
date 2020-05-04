module Elm exposing (main)

import Browser
import Html exposing (Html)
import Html.Attributes exposing (href, for, class, placeholder)

import Bootstrap.CDN as CDN
import Bootstrap.Form as Form
import Bootstrap.Form.Textarea as Textarea
import Bootstrap.Button as Button
import Bootstrap.Popover as Popover
import Bootstrap.Alert as Alert
import Bootstrap.Grid as Grid
import Bootstrap.Grid.Row as Row
import Bootstrap.Grid.Col as Col
import Bootstrap.Text as Text


type OperationType = Contigs | Peptides

type alias Model =
    { optype : Maybe OperationType
    , facontent : String
    , helpPopoverState : Popover.State
    }

type Msg = NoMsg
        | SelectOp OperationType
        | UpdateFacontent String
        | HelpPopover Popover.State
        | SetExample
        | SubmitData


main : Program () Model Msg
main =
  Browser.document
    { init = init
    , view = view
    , update = update
    , subscriptions = \_ -> Sub.none
    }

init : () -> (Model, Cmd Msg)
init () = (
        { optype = Nothing
        , facontent = ""
        , helpPopoverState = Popover.initialState
        }, Cmd.none)

update : Msg -> Model -> (Model, Cmd Msg)
update m model = case m of
    NoMsg -> (model, Cmd.none)
    SelectOp p ->
        -- Iff the example input is selected, switch it
        if model.optype == Just Contigs && model.facontent == contigExampleData && p == Peptides then
            ({ model | optype = Just Peptides, facontent = peptidesExampleData }, Cmd.none)
        else if model.optype == Just Peptides && model.facontent == peptidesExampleData && p == Contigs then
            ({ model | optype = Just Contigs, facontent = contigExampleData }, Cmd.none)
        else ({ model | optype = Just p }, Cmd.none)
    UpdateFacontent c -> ({model | facontent = c}, Cmd.none)
    HelpPopover state -> ({ model | helpPopoverState = state}, Cmd.none )
    SetExample ->
        let
            nc = case model.optype of
                Nothing -> "?" -- should never happen
                Just Contigs -> contigExampleData
                Just Peptides -> peptidesExampleData
        in ({model | facontent = nc}, Cmd.none)
    SubmitData -> ( model, Cmd.none ) -- TODO


validateFasta : OperationType -> String -> Maybe String
validateFasta p fa =
    if String.length fa > 50200 -- 200 is a little margin
        then Just "Input is too large! Only the first 50,000 characters will be analyzed. Run the tool locally to remove any limitations"
        else if nrSeqs fa > 1004 -- 4 is a little margin
            then Just "Too many sequences. Only the first 1,000 sequences will will be analyzed. Run the tool locally to remove any limitations"
            else case p of
                Peptides -> Nothing
                Contigs -> Nothing -- validate that it's not AAs

nrSeqs : String -> Int
nrSeqs fa = String.filter (\c -> c == '>') fa |> String.length

view : Model -> Browser.Document Msg
view model =
    { title = "AMP Prediction"
    , body = [ CDN.stylesheet
             , CDN.fontAwesome
             , layout model ]
    }


layout model =
    Grid.container []
    [ Grid.simpleRow
        [Grid.col []
            [ header
            , Html.hr [] []
            , intro
            , Html.hr [] []
            , viewModel model
            , Html.hr [] []
            , outro
            , Html.hr [] []
            , footer
            ] ]

    ]

header : Html Msg
header = Grid.simpleRow
         [Grid.col [] [Html.text "Macrel"]
         ,Grid.col [] [Html.text "Prediction"]
         ,Grid.col [] [Html.text "About"]]

intro : Html Msg
intro =
    Alert.simpleInfo [] [ Html.p [] [Html.text "If you use macrel in your published work, please cite:"]
                         , Html.blockquote [] [Html.p []
                                    [Html.em [] [Html.text """
                                    MACREL: antimicrobial peptide screening in genomes and metagenomes. Celio Dias Santos-Junior, Shaojun
                                    Pan, Xing-Ming Zhao, Luis Pedro Coelho. bioRxiv 2019-2020; DOI:"""
                                    ,Html.a [href "https://doi.org/10.1101/2019.12.17.880385"] [Html.text "10.1101/2019.12.17.880385"]]]
                        , Html.p [] [Html.text "The preprint also details the macrel algorithms and presents benchmarking results"]]]

outro : Html Msg
outro = Html.div []
    [Html.p [] [Html.text """Macrel uses machine learning to select peptides
    with high probability of being an AMP. Macrel is optimized for higher
    specificity (low rate of false positives)."""]
    ,Html.p [] [Html.text """Macrel will also classify AMPs into hemolytic and
    non-hemolytic peptides."""]
    ]

footer : Html Msg
footer = Html.text "Copyright 2019-2020 Macrel authors"

viewModel : Model -> Html Msg
viewModel model =
    let
        buttonStyle who active = case active of
            Nothing -> [ Button.primary , Button.onClick (SelectOp who)]
            Just p -> if who == p
                        then [ Button.info, Button.onClick (SelectOp who) ]
                        else [ Button.outlineSecondary , Button.onClick (SelectOp who)]
        placeholderText = case model.optype of
            Nothing -> "Select input type above..."
            Just Contigs -> ">ContigID\nAATACTACTATCTCTCTCTACTATCTACATCATCA...\n"
            Just Peptides -> ">PeptideID\nMEPEPAGAD....\n"
        faerror = case model.optype of
            Nothing -> Nothing
            Just p -> validateFasta p model.facontent

    in Grid.simpleRow
        [ Grid.col [] <|
            [Html.h2 [] [Html.text "Online AMP prediction"]
            , Html.p [] [Html.strong [] [Html.text "Step 1."]
                        ,Html.text " Select mode:"]
            , Grid.simpleRow
                    [ Grid.col [] [Button.button (buttonStyle Contigs model.optype) [ Html.text "Predict from contigs (DNA sequences)" ]]
                    , Grid.col [] [Button.button (buttonStyle Peptides model.optype) [ Html.text "Predict from peptides (amino acid sequences)" ]]]
            , Html.p [] [Html.text
                "(The command line tool also supports prediction from short-reads, but this is not available on the webserver)."]
            , case faerror of
                    Nothing -> Html.text ""
                    Just err -> Alert.simpleWarning [] [ Html.text err]
            , case model.optype of
                Nothing -> Html.text ""
                Just p ->
                    Form.group []
                        [ Html.label [ for "fasta" ]
                                    [ Html.strong [] [Html.text "Step 2."]
                                    , Html.text <| if p == Contigs
                                                        then " Input DNA FASTA "
                                                        else " Input Peptides FASTA "
                                    , Popover.config
                                        ( Button.button
                                            [ Button.small
                                            , Button.primary
                                            , Button.attrs <|
                                                Popover.onHover model.helpPopoverState HelpPopover
                                            , Button.attrs <|
                                                Popover.onClick model.helpPopoverState HelpPopover
                                            ]
                                            [ Html.span [class "fa fa-question-circle"] [] ]
                                        )
                                        |> Popover.right
                                        |> Popover.titleH4 [] [ Html.text "FASTA format" ]
                                        |> Popover.content []
                                            [ Html.text (case model.optype of
                                                        Nothing -> ""
                                                        Just Contigs -> "Please provide nucleotides (or change to peptides mode above)."
                                                        Just Peptides -> """
Peptides submitted to the Macrel prediction should consist of 20 canonical
amino acids and their length should range from 10 to 100 amino acids. Please
avoid contigs containing non-canonical bases, such as N, R or Y.""") ]
                                        |> Popover.view model.helpPopoverState
                                    ]
                            , Textarea.textarea
                                <| [ Textarea.id "fasta"
                                , Textarea.rows 10
                                , Textarea.onInput UpdateFacontent
                                , Textarea.attrs [placeholder placeholderText]
                                , Textarea.value model.facontent
                                ] ++ (case faerror of
                                        Nothing -> []
                                        Just _ -> [Textarea.danger])
                        , Grid.row [Row.rightXl]
                            [ Grid.col [] [Html.text ""]
                            , Grid.col [Col.textAlign Text.alignXsRight]
                                [Button.button [ Button.small, Button.outlineSecondary, Button.onClick SetExample ] [Html.text "Example"]]]
                        , Button.button [ Button.primary, Button.onClick SubmitData] [ Html.text "Submit" ]
                        ] ] ]


contigExampleData : String
contigExampleData = """>scaffold2530_2_MH0058
CTTCTGATCTTTACGCAGCATTGTGTGTTTCCACCTTTCAAAAAATTCTCCGTGAACTGC
GCCCTGGGAGTGGTGAAATCCTCCGCGGAACGAAGTCCCGGAATTGCGCACAAATTCACG
TGCTGAACAATTTTACCATAGGAATGTGCGGTTGTAAAGAGAAAAATGCAAAAAATTCCT
TATTTTTATAAAAGGAGCGGGGAAAAGAGGCGGAAAATATTTTTTTGAAAGGGGATTGAC
AGAGAGAAACGGCCGTGTTATCCTAACTGTACTAACACACATAGTACAGTTGGTACAGTT
CGGAGGAACGTTATGAAGGTCATCAAGAAGGTAGTAGCCGCCCTGATGGTGCTGGGAGCA
CTGGCGGCGCTGACGGTAGGCGTGGTTTTGAAGCCGGGCCGGAAAGGAGACGAAACATGA
TGCTGTTTGGTTTTGCGGGGATCGCCGCCATCGTGGGTCTGATTTTTGCCGCTGTTGTTC
TGGTGTCCGTGGCCTTGCAGCCCTGAGAACGGGGCAGATGCAATGAGTACGCTGTTTTTG
CTTGGTATCGCAGGCGCGGTACTGCTGGTTATTTTGCTGACAGCGGTGATCCTGCACCGC
TGATCGAACATTCCTCAGAAAGGAGAGGCACACGTTCTGACATTGAATTACCGGGATTCC
CGTCCCATTTATGAACAGATCAAGGACGGCCTGCGGCGGATGATCGTCACCGGGGCC
>scaffold75334_1_MH0058
ACAGCTTCCTCACCATCAACAGCCACTGCTACGATACCGGCAGGAACAGAGATTGTAGCG
TTATCGGAAGTAAGAACGGTCTCAGCGATTACCTTACCCAAGATATTCGTGATAACTACA
GACTTACCAGCTGCGCCTTGAACGGTTACGGTACCGTTACCGGCTACTACAGAGATACCT
TCTACAGCATCGATTGTTTCGTTATCTGTTGCGATATCATCACCTTGCTCCACATTAAAG
ATCAAAGCATCGTCACCACCTGTATTCATCTCATCGAAATTAGAAGAACGATCATCAGAC
AGAACTAAACAACCATTCTGCATTCTCAACCAAGCCGCATATGTAGGAGCGATGTCACCC
ATAGCAGAATTATTGTAACCAGGAACATGTTTACTCTTCATAGACTCGAACAAGAATGCT
CTATCAGCTTCTACCTCATTAGCGGCAACTTCCTTGTTTACGAAACGCATAGACCAAGTC
ACATACTTATGGTTATCACCAGATAAGATATATTTGTGTGAAACAAATTCATCTTCAGAT
ACTCCAGCCTTCTTAGCCGCAGTCCAAGCATCCTCCTCAGCCTTGTTCAAATCAGCGAAG
TTGATCTTCTCGTTCGGTAAGTTCTTGAACTCATCTCTCAAGATATACAAAGTATCAGCT
ACACGGATAGCTTCCACGAAACCTGCACGATCATATTTCTTCCACATGTAGTCAGCATCT
TGCTCCCCATTTGATAATGTTGACCAACCACATTTATTAGCGAAATCATGGAAATTGATC
AAATATTTACCTCTTTCAAAGCCCGGAACAGCCGGTTTTGCGTGAACGCAATGCCATTTA
TCTGTCGGTTTACCTTCTGCGGTAAAGTGATGGTCATCTTCCGTACAAGGTACGCCTTCA
ACTCCTTCGAAATCGTTACGATCAATAGAGATCAAATATTGAGGCTTGATATTACCAGAA
CCACGGTTTACCCATGCGGTATCAACGATGAATGACAAACCATCCTCAGTCTTATCCGGA
GTATAAATACCTAAGAAGTCAATACCTTCTTTCATGAAGTTCTTATTATTCTCAACCTGC
AAGTATTCTTTACGGTATTTCTCGATAAAACGTAAAGTATCGGCCTTGTCACCTTCATTA
CCTTCAAGTTCAAGAGAACTGAAACGACGGTAAAGCGGAGTGTTGTCCGGTTCGATAGCG
AAAGCAGAAGTACGAGTCTCGCCTAATACTTGATTCTTCAATGTAGCGGCACCGTCATAA
TCGGATACTCCAGCTTTCTGATAACCTAATACATATTTTGTAGCATTAGAATAATTATCA
TATGCAGCATTAACAATCGCATAATAGTGCTTTCCACTGATATAGTTGTTTTCTTTGAAG
AATGCATAATAAGTCGTATAAGATTGCTGATGACCAGTAGCACTAGAACCTACATTAAAC
TTATCTTCCTTATTAACTTTCCATTCATTAAGTTCACCTTTACCATCTTTATAAGAAACC
TCATAAGCAGTTCTCACTAGAGTCTTCAAGCCTTTAATACGATTCTTAGCCACATCATCG
CTAACCTTATAACCGTAAGGGATATCTACAGTCGCAAATATCTTTGAGTTATTATTCCAA
CGCTCAGTCAATGTATCAATTACGAAAGCATCTTTTCCATCCAACACAGTCAACGTAGAG
TCTGTAGATTTAGCGATGTATTTATCCGTAGCATAAGGATGCCAATAGTTAAACGTATAT
TTGTTAACCAACAAAGAATCTTTCTCCAACTTTCTGTAACCTAAATACGGATCAGAAATA
GAAGCTGCCGGAACTTTCTCGAAGAACAAAGAATCTATACCGACAGAGCCATACTCCGCT
GTAAGTTTATCTCCATATACAAACATATAAGAACCACCCTCTGCCTTTCTTAACTGAACA
GTAGAATAAACAGCAGTCTTATATGACTCCATCTGGCCATCACCATCAAAATCCGCATTT
ATAACTTTATCTGCGAATTCACGGTTAGAGATAGTAACCGGAGACACAGCCTTCACTTTA
TCATTACTTGTATTTTTCTTCAATACAACCCATTGATAAGCTGGCATATGAGCAACACTC
TGCTCATCCTCATTCACGGTTGTCCAACGGATAGTACCATTCTCATAGATAGGAGAAGCT
AAATACTGTCCTTGTTTATTCTTGATGAAATAAACACCATCATCAACGGAAGTCTTGGAA
GAACCCTTTTCCTCACATCCGCTAAGACCCAAAGAAATATGAGTGTTCTGGTCCTCCGGA
TCAACTGTCAAAATACAAGTTTCATCTTTAATCAAGTCTTGCAAAGAGACACGCCAGAAC
TCATGCTCGTTAACTGGAGTAGCCGTAACGTGCGTATTCCAATACTTGCCATCATTATTC
CAAGTAACTTTTTTTACATAGATGTATAAGCTATCGCCACTCGGAGAATAGATGAATTTG
AACTGATGTTGAGCACGTAATTCAGAGCTAATAGAACTGTAAGCATCAGCTTTATCCTTT
GCTTTCTCTGTCCAACCGTAAGCTAAGAACTTTGTACCTGTCTCATTCGTATAAGCCGTA
TCAACTTTCAAATAAGCGTCTTTATCCTTTTGCTTGATGAAAAGCCATTTGTTATTATCT
TTATCCTCTGCGATAAACTTCTGCTCATTGAACGGATTCTTTAAGCTAGTACCTTTCACG
TCCGGAGTGAAAGTCAATTGTACTCCAGCCTTGTTATCTTGAAGAAGACCCAACTTCGTA
TTTACTTGATCTTCGCTCAATGCGATTTGTGCTGCGCCAACCAAAGCAAAGTTAGTCACA
TCAGCGGCGTTTGAAGTGGCATCTATTTGATTAGCACCCCACTTGGCAACCTTAACAGCT
CCGGTAGTAGGATCCGCTTTCAAACCAACAATTGAGTCCGTTGAGAAATAAGAATACAAA
GGTCTCTTTTCTTCCAAATCCTTGTATGAGCGAGAGAATGCCCAACCTGCGATCTCACCA
CCTAAAACAGGTTTCCAAGCACCATTCGTACCATTCTCAGTCTTCTCATGGCCAGCCATT
GTCAAATCCAACATAGTTCCAGACAATTTATTCTGGAAGTCATAAATAGGAGCTTGACCT
TGATTATAATTAGAGACAGAAACACACCACAATGTAGCCTCTAAGTTATTTTTAGCTTCG
CTTGCACTATAAACACGTAGTTCATAATCACCTGTTCCACGATTTAGCTCCATCGCAAGA
TAAGCAGGAGTAGTGCCATCCATGACCTTCAACTGATAAAGACCAGAATTAGCTCCCTCC
TTCAAACCGCCTAAACGCCATTCTTGACCCAATACAGTTTCTGGGTCTACTCGACTTGGT
GTAGTCTGTGCATTAACAGACATAACACTTAACAGTGCCATACCTGCCAAAAGAGTAGAA
AACT"""

peptidesExampleData : String
peptidesExampleData = """>AP00002|AMP
YVPLPNVPQPGRRPFPTFPGQGPFNPKIKWPQGY
>AP00007|AMP
GNNRPVYIPQPRPPHPRL
>P20491|NAMP
MISAVILFLLLLVEQAAALGEPQLCYILDAVLFLYGIVLTLLYCRLKIQVRKAAIASREKADAVYTGLNTRSQETYETLKHEKPPQ
>P35109|NAMP
MVDDPNKVWPTGLTIAESEELHKHVIDGSRIFVAIAIVAHFLAYVYSPWLH
>P19962|NAMP
DSDSAQNLIG"""

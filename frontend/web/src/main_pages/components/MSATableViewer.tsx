import * as React from "react"; 

export interface MSATAbleViewerProps { alignments: any }

export default class MSATableViewer extends React.Component<MSATAbleViewerProps, {}> {

    xrefcaptions: any;

    constructor(props: MSATAbleViewerProps) {
        super(props);

        this.xrefcaptions = {
            'GO': 'Gene Ontology (GO)',
            'Pfam': 'Protein Families',
            'Interpro': 'Interpro'
        };
    }

    componentDidMount()
    {
    }

    render() {

        var allRows = []
        var allXrefs = []

        allXrefs.push(<tr key={0}>
        <th colSpan={2} >Name</th>
        <th>{this.xrefcaptions.GO}</th>
        <th>{this.xrefcaptions.Pfam}</th>
        <th>{this.xrefcaptions.Interpro}</th>
        </tr>)

        for (var i = 0; i < this.props.alignments.align.length; ++i)
        {
            var rowAlign = this.props.alignments.align[i];

            allRows.push(
            <tr key={i}>
                <th style={{ position:"absolute", left:0, width: "100px", verticalAlign: "top", borderTop: "1px solid #ccc", padding: "10px"}}>{rowAlign.id}</th>
                <td style={{ background:"yellow", verticalAlign: "top", borderTop: "1px solid #ccc", padding: "10px"}}><pre>{rowAlign.alignment}</pre></td>
            </tr>
            );

            var xKeys = Object.keys(rowAlign.xrefs);

            var xrefvals = [];

            for (var l=0; l < xKeys.length; ++l)
            {
                var xkey = xKeys[l];
                var xKeyVals = rowAlign.xrefs[xkey];

                for (var xv=0; xv < rowAlign.xrefs[xkey].length; ++xv)
                {
                    xrefvals.push( rowAlign.xrefs[xkey][xv] )
                }
            }

            allXrefs.push(<tr key={allXrefs.length}>
                <td>{rowAlign.org}</td>
                <td>{rowAlign.seqid}</td>

                <td>{rowAlign.xrefs['GO'].map((x, i) => <a key={i} style={{marginLeft: "10px"}} href={"http://amigo.geneontology.org/amigo/term/"+x}>{x}</a>)}</td>
                <td>{rowAlign.xrefs['Pfam'].map((x, i) => <a key={i} style={{marginLeft: "10px"}} href={"https://pfam.xfam.org/family/"+x}>{x}</a>)}</td>
                <td>{rowAlign.xrefs['Interpro'].map((x, i) => <a key={i} style={{marginLeft: "10px"}} href={"https://www.ebi.ac.uk/interpro/entry/"+x}>{x}</a>)}</td>

            </tr>);

        }

        return (
            <div>
            <h1>{this.props.alignments.id}</h1>

            <div style={{position:"relative"}}>
                <div style={{overflowX:"scroll", overflowY:"visible", width:"90%", marginLeft:"110px"}}>
                    <table style={{ tableLayout: "auto", width: "100%"}}><tbody>
                        {allRows}
                        </tbody>
                    </table>
                </div>
        </div>
        <div>
            <h3>Cross-References</h3>
        <table style={{ tableLayout: "auto", width: "100%"}}><tbody>
                        {allXrefs}
                        </tbody>
                    </table>
            </div>
        </div>);
    }
}
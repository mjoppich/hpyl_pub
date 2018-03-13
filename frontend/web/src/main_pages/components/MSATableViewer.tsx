import * as React from "react"; 
import FlatButton  from 'material-ui/FlatButton';
import DownloadArchive from 'material-ui/svg-icons/content/archive';

export interface DownloadButtonProps { filename: string, getDownloadContent: any }


class DownloadButton extends React.Component<DownloadButtonProps, {}> {

    makeAlignmentDownload()
    {
        var fileContents = this.props.getDownloadContent();

        var filetype = "text/plain";

        var a = document.createElement("a");
        var dataURI = "data:" + filetype +
            ";base64," + btoa(fileContents);
        a.href = dataURI;
        a['download'] = this.props.filename;
        var e = document.createEvent("MouseEvents");
        // Use of deprecated function to satisfy TypeScript.
        e.initMouseEvent("click", true, false,
            document.defaultView, 0, 0, 0, 0, 0,
            false, false, false, false, 0, null);
        a.dispatchEvent(e);

        a.parentNode.removeChild(a);
    }

    render()
    {
        return <FlatButton 
        label=""
        labelPosition="before"
        primary={true}
        icon={<DownloadArchive />}
        onClick={() => this.makeAlignmentDownload()}
        />
        
    }

}

export interface MSATAbleViewerProps { alignments: any }

export default class MSATableViewer extends React.Component<MSATAbleViewerProps, {}> {

    xrefcaptions: any;
    goNS2Name: any;

    constructor(props: MSATAbleViewerProps) {
        super(props);

        this.xrefcaptions = {
            'GO': 'Gene Ontology (GO)',
            'TMGO': 'Textmining GO',

            'Pfam': 'Protein Families',
            'Interpro': 'Interpro'
        };

        this.goNS2Name = {
            'biological_process': 'Biological Process',
            'molecular_function': 'Molecular Funciton',
            'cellular_component': 'Cellular Component'
        }
    }

    getGONS2Name( namespaceID, defaultVal )
    {
        if (this.goNS2Name[namespaceID] !== undefined)
        {
            return this.goNS2Name[namespaceID];
        }

        return defaultVal;
        
    }


    componentDidMount()
    {
    }

    makeGOLinks( rowAlign, idCat )
    {
        var allGOClasses = Object.keys(rowAlign.xrefs[idCat])

        console.log("all GO Classes");
        console.log(allGOClasses);
        console.log(allGOClasses.map((namespaceID, nsi) => rowAlign.xrefs[idCat][namespaceID]));

        var self = this;


        var elemsByNamespace = allGOClasses.map((namespaceID, nsi) => 
            <div key={nsi}>
                <h3>{self.getGONS2Name(namespaceID, namespaceID)}</h3>
                {rowAlign.xrefs[idCat][namespaceID].map(
                    (x, i) => 
                                <div key={i}>
                                    <a href={"http://amigo.geneontology.org/amigo/term/"+x}>{x}</a><br/>
                                    <span>{rowAlign.xrefs["GOTERMS"][x]}</span><br/>
                                </div>
                )}
            </div>
        );

        return <div>{elemsByNamespace}</div>
    }

    printMSA()
    {
        var outStr = "";
        for (var i = 0; i < this.props.alignments.align.length; ++i)
        {
            var rowAlign = this.props.alignments.align[i];

            var rowID = rowAlign.entryID;
            var alignSeq = rowAlign.alignment;

            outStr += ">" + rowID + "\n" + alignSeq + "\n";
        }

        return outStr;
    }

    printSeq( seqid )
    {
        var outStr = "";
        for (var i = 0; i < this.props.alignments.align.length; ++i)
        {
            var rowAlign = this.props.alignments.align[i];

            var rowID = rowAlign.entryID;
            var alignSeq = rowAlign.alignment;

            if (rowID == seqid)
            {
                outStr += ">" + rowID + "\n" + alignSeq + "\n";
            }
        }

        return outStr;
    }

    makeProtIDInfo(rowAlign)
    {

        let thisRowID = rowAlign.entryID;

        var recID = rowAlign.recordID;
        var recStart = rowAlign.genomicStart;
        var recEnd = rowAlign.genomicEnd;
        var recStrand = "?";
        
        if (rowAlign.genomicStrand == -1)
        {
            recStrand = "-";
        } else if (rowAlign.genomicStrand == 1)
        {
            recStrand = "+"
        }

        var posInfo = recID + " " + recStart + "-" + recEnd + ":" + recStrand;

        return <div><p>{rowAlign.entryID}</p><DownloadButton filename={thisRowID + ".fa"} getDownloadContent={() => this.printSeq(thisRowID)}/><p>{posInfo}</p></div>
    }

    render() {

        var allRows = []
        var allXrefs = []

        allXrefs.push(<tr key={0} style={{textAlign: 'left'}}>
        <th style={{width: "100px"}}>Organism</th>
        <th style={{width: "150px"}}>Gene/Protein Name</th>
        <th style={{width: "100px"}}>UniProt</th>

        <th>{this.xrefcaptions.GO}</th>
        <th>{this.xrefcaptions.TMGO}</th>

        <th>{this.xrefcaptions.Pfam}</th>
        <th>{this.xrefcaptions.Interpro}</th>
        </tr>)

        for (var i = 0; i < this.props.alignments.align.length; ++i)
        {
            var rowAlign = this.props.alignments.align[i];

            allRows.push(
            <tr key={i}>
                <th style={{ position:"absolute", left:0, width: "100px", verticalAlign: "top", borderTop: "1px solid #ccc"}}>{rowAlign.entryID}</th>
                <td style={{ background:"yellow", verticalAlign: "top", borderTop: "1px solid #ccc"}}><pre style={{margin: "0"}}>{rowAlign.alignment}</pre></td>
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
                <td style={{verticalAlign: 'top'}} ><p>{rowAlign.organismName}</p><p>{rowAlign.organismID}</p></td>
                <td style={{verticalAlign: 'top'}}>{this.makeProtIDInfo(rowAlign)}</td>
                <td style={{verticalAlign: 'top'}}>{rowAlign.xrefs['Uniprot'].map((x, i) => <div key={i}><a href={"http://www.uniprot.org/uniprot/"+x}>{x}</a><br/></div>)}</td>

                <td style={{verticalAlign: 'top'}}>{this.makeGOLinks(rowAlign, 'GO')}</td>
                <td style={{verticalAlign: 'top'}}>{this.makeGOLinks(rowAlign, 'TMGO')}</td>

                <td style={{verticalAlign: 'top'}}>{rowAlign.xrefs['Pfam'].map((x, i) => <div key={i}><a href={"https://pfam.xfam.org/family/"+x}>{x}</a><br/></div>)}</td>
                <td style={{verticalAlign: 'top'}}>{rowAlign.xrefs['Interpro'].map((x, i) => <div key={i}><a href={"https://www.ebi.ac.uk/interpro/entry/"+x}>{x}</a><br/></div>)}</td>

            </tr>);

        }

        return (
            <div>
            <h1>{this.props.alignments.id}
            <DownloadButton filename={this.props.alignments.id + ".fa"} getDownloadContent={() => this.printMSA()}/>

                </h1>
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
        <table style={{ tableLayout: "fixed", width: "100%"}}><tbody>
                        {allXrefs}
                        </tbody>
                    </table>
            </div>
        </div>);
    }
}
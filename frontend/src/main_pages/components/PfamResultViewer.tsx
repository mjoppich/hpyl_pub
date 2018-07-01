import * as React from 'react'; 

import ChipInput from 'material-ui-chip-input'
import AutoComplete from 'material-ui/AutoComplete'
import axios from 'axios';
import config from '../config';

import SelectField from 'material-ui/SelectField';
import MenuItem from 'material-ui/MenuItem';

export interface PfamResultViewerProps { pfamres: any};
export interface PfamResultViewerState { selectedElement: any, resultPart: any};

export default class PfamResultViewer extends React.Component<PfamResultViewerProps, PfamResultViewerState>{

    constructor(props)
    {
        super(props);
    }

    componentWillMount()
    {
        var self=this;
        this.setState({selectedElement: undefined, resultPart: ""});

    }

    componentDidMount(){

    }

    handleChange(event, index, value)
    {
        var newValue = value;
        this.setState({selectedElement: newValue})

        var table = this.makeResult(newValue);
        this.setState({resultPart: table});
    }

    makeResult(allowedGeneIDs: Array<any>)
    {
        console.log("make result");
        console.log(allowedGeneIDs);

        var errorMessages = [];
        var structStuff = [];


        for (var u = 0; u < this.props.pfamres.length; ++u)
        {
            var resultElem = this.props.pfamres[u];

            if (allowedGeneIDs.indexOf(resultElem['GENE_ID']) == -1)
            {
                continue;
            }

            /*
                    pfamRes= {
                        'GENE_ID': seqID,
                        'GENE_START': alignStart,
                        'GENE_END': alignEnd,
                        'ENV_START': envStart,
                        'ENV_END': envEnd,
                        'PFAMID': pfamID,
                        'PFAMNAME': pfamName,
                        'PFAM_TYPE': pfamType,
                        'HMM_START': pfamHMMStart,
                        'HMM_END': pfamHMMEnd,
                        'HMM_LENGTH': pfamHMMLength,
                        'PFAM_BIT_SCORE': pfamBitScore,
                        'PFAM_EVALUE': pfamEvalue,
                        'PFAM_SIGNIFICANT': pfamSignificant
                    }

                    if pfamClan:
                        pfamRes['PFAM_CLAN'] = pfamClan

                    if hasActiveSites:
                        pfamRes['PFAM_ACTIVE_SITES'] = pfamActiveSites
            
            */

            
            var evalue = parseFloat((resultElem['PFAM_EVALUE']).toFixed(5));
            var bitscore = parseFloat((resultElem['PFAM_BIT_SCORE']).toFixed(5));

            var geneID = resultElem['GENE_ID'];
            var geneRegion = resultElem['GENE_START'] + "-" + resultElem['GENE_END'];
            var envRegion = resultElem['ENV_START'] + "-" + resultElem['ENV_END'];
            var hmmRegion = resultElem['HMM_START'] + "-" + resultElem['HMM_END'] + " ("+resultElem['HMM_LENGTH']+" AAs)";

            var pfamName = resultElem['PFAMNAME'];
            var pfamID = resultElem['PFAMID'];

            var activeSites = []
            if ('PFAM_ACTIVE_SITES' in resultElem)
            {
                activeSites = resultElem['PFAM_ACTIVE_SITES'];
            }

            var clanName;
            if ('PFAM_CLAN' in resultElem)
            {
                var cn = resultElem['PFAM_CLAN'];
                clanName = <a target="_blank" href={"https://pfam.xfam.org/clan/"+cn}>{cn}</a>;
            } else {
                clanName = "";
            }

            var rowStyle = {};
            if (structStuff.length % 2 === 0)
            {
                rowStyle = {
                    backgroundColor: config.table_even_bg
                }
            }

            structStuff.push(
                <tr key={structStuff.length} style={rowStyle}>
                    <td style={{verticalAlign: 'top'}}>{geneID}</td>
                    <td style={{verticalAlign: 'top'}}>{geneRegion}<br/>{envRegion}</td>
                    <td style={{verticalAlign: 'top'}}>{hmmRegion}<br/>{activeSites.join(", ")}</td>
                    <td style={{verticalAlign: 'top'}}>{bitscore}<br/>{evalue}</td>

                    <td style={{verticalAlign: 'top'}}>
                        <a target="_blank" href={"https://swissmodel.expasy.org/templates/"+pfamID}>{pfamID}</a><br/>
                        {pfamName}<br/>
                        {clanName}
                        </td>
                </tr>
            );
        }

        var tableElems = <div></div>;

        if (structStuff.length > 0)
        {
            tableElems = <table style={{ tableLayout: "unset", width: "100%"}} cellSpacing={0} cellPadding={0} >
                            <tbody>
                                <tr style={{textAlign: 'left'}}>
                                    <th>Gene ID</th>
                                    <th>Gene Region<br/>Envelope Region</th>
                                    <th>HMM Region <br/>Predicted Active Sites</th>
                                    <th>Pfam Bit Score <br/>Pfam E-Value</th>
                                    <th>Pfam Name<br/>Pfam ID<br/>CLAN</th>
                                </tr>
                                {structStuff}
                            </tbody>
                        </table>;
        }

        return (
            <div>
            {tableElems}
            {errorMessages}
            </div>
        );


    }
  
    render(){

        console.log("PfamResultViewer")
        console.log(this.props.pfamres);

        var menuitems = [];

        var allPfamSeqs = [];


        for (var i = 0; i < this.props.pfamres.length; ++i)
        {
            
            var elem = this.props.pfamres[i];
            var seqid = elem['GENE_ID'];

            if (allPfamSeqs.indexOf(seqid) == -1)
            {
                allPfamSeqs.push(seqid);

                var elemitem = <MenuItem key={i} value={seqid} primaryText={seqid} />
                menuitems.push(elemitem);
            }
        }

        var menuPart =  <SelectField hintText="Select one or multiple proteins for Pfam lookup." multiple={true} value={this.state.selectedElement} onChange={this.handleChange.bind(this)}>
                            {menuitems}
                        </SelectField>

        return (<div>
            
            <h3>Select the protein:</h3>
            {menuPart}
            {this.state.resultPart}

            </div>);
    }
}